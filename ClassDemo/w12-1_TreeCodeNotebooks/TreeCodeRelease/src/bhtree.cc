/*
  Time-stamp: <bhtree.cc on Monday, 15 November, 2021 at 23:22:20 MST (pinto)>

  Simple Barnes-Hut tree implementation based on the octree in octree.cc.

  Several different traversal algorithms are included at the bottom of this file

*/

#ifndef  __BHTREE_INCLUDE__
#define  __BHTREE_INCLUDE__

#include "avxdirect.cc"
#include "octree.cc"

// execute parallel "code" if getStats is true
#define STATS(code) if (getStats) { int g = omp_get_thread_num(); code }

template <typename T>
class BHtree : public Octree<T> {
public:

    using vec = SmallVec<T,3>;

    // Barnes-Hut node data
    struct bhdata {
        vec com;             // centre of mass
        T mass;              // mass of node
        T Bmax;              // maximum distance of particle from com
        double B2, B3, B4;   // parameters in error bound
        double extent;
    };

    // functions from octree used herein:
    using Octree<T>::buildTree;
    using Octree<T>::isLeaf;
    using Octree<T>::nodeSize;
    using Octree<T>::begin;
    using Octree<T>::end;
    using Octree<T>::isIn;
    using Octree<T>::pindex;
    using Octree<T>::tree;
    using Octree<T>::ps;
    using Octree<T>::np;
    using Octree<T>::nNodes;
    using Octree<T>::ROOT;
    using Octree<T>::myprint;
    using Octree<T>::putParticlesInOriginalOrder;

    BHtree( ParticleSet<T> &pset, int maxLeafSize, T epsSmooth, int maxSources);
    ~BHtree();

    // tree-building
    void makeTree();
    void propagateCOM(int node);
    void checkCOM(int node);
    void propagateBmax(int node);
    void propagateExtent(int node);

    // traversals
    void BHsubsets(T theta, int maxILlength, vec delta);
    void recursiveSubsetPartition(T theta, int maxILlength, vec delta);
    void acc1(int p, T theta, vec delta);
    void acc1_noAVX(int p, T theta, vec delta);
    void accAll(T theta, vec delta);
    void accAll_noAVX(T theta, vec delta);

    // internals for traversals
    void direct_noAVX(int sink, int srcBeg, int srcEnd, vec delta);
    void direct_noAVX(int sink, vec source, T srcmass, vec delta);

    void bhSingle(int p, T theta, int node, vec delta);
    void bhSingle_noAVX(int p, T theta, int node, vec delta);

    void bhLeaf(int sinknode, vec delta, T theta);
    void recursiveSubsetPartitionInternal(int node, T theta, int maxILlength, vec delta, int level);
    void bhTraverse(T BmaxA, vec posA, int nodeB, T theta, int maxmax, bool flag, int thread);

    void BHsubsetsInternal(int sinknode, T theta, int maxILlength, vec delta, int level);

    bool CellAContainsCellB(int A, int B);
    bool CellIntersectsBall(int Cell, vec &centre, T radius);

    bool BHThetaMAC( int p, T theta, int node );

    void resetStats();
    void reportStats();

    bhdata *BH;
    T eps2;
    AVXDirect<T> *avx;
    int nprocs;

    bool getStats = false;
    uint64 *nleaf, *nnode, *nbar;
    uint64 *ndirect, *illn, *illl;
};

template <typename T>
BHtree<T>::BHtree( ParticleSet<T> &pset, int maxLeafSize, T epsSmooth, int maxsrclen)
                   : Octree<T>(pset, maxLeafSize),  eps2( epsSmooth * epsSmooth) {
    BH = NULL;
    nprocs = omp_get_num_procs();
    // allocate nprocs copies of avxDirect
    alloc(avx, nprocs);
    for(int g=0; g<nprocs; g++) avx[g].initAVX(maxsrclen);

}

template <typename T>
BHtree<T>::~BHtree() {
    if( BH != NULL ) delete[] BH;
    delete[] avx;
    STATS(
        delete[] nleaf;   delete[] nnode;
        delete[] nbar;    delete[] ndirect;
        delete[] illn;    delete[] illl;
    );
}



template <typename T>
void BHtree<T>::resetStats() {
    STATS(
        for(int g=0; g<nprocs; g++) {
            nleaf[g] = 0;   nnode[g] = 0;
            nbar[g]  = 0; ndirect[g] = 0;
            illn[g]  = 0;    illl[g] = 0;
        }
        );
}

template <typename T>
void BHtree<T>::reportStats() {

    STATS(
            // sum partial sums of statistics over threads
            for(int i=1; i<nprocs; i++) {
                nleaf[0]+= nleaf[i];    nnode[0]+= nnode[i];
                illl[0]+=illl[i];       illn[0]+=illn[i];
                ndirect[0]+=ndirect[i]; nbar[0]+=nbar[i];
            }
            double avgnbar = (double)nbar[0]/(double)(nnode[0]+nleaf[0]);

            myprint("           number of leaf cells as sinks: %d\n", nleaf[0]);
            myprint("           number of node cells as sinks: %d\n", nnode[0]);
            myprint("                    averge sink set size: %7.2e\n", avgnbar);

            double illlbar = 0;
            if(nleaf[0]>0) illlbar = (double)illl[0]/(double)nleaf[0];
            myprint("<interaction list length for leaf sinks>: %7.2e\n", illlbar);

            double illnbar = 0;
            if(nnode[0]>0) illnbar = (double)illn[0]/(double)nnode[0];
            myprint("<interaction list length for node sinks>: %7.2e\n", illnbar);

            myprint("                 total number of directs: %7.2e\n", ndirect[0]);
            myprint("                         fraction vs N^2: %7.2e\n", (double)ndirect[0]/((double)np*(double)np));
            myprint("\n");
          );

        resetStats();
}

template <typename T>
void BHtree<T>::makeTree( ) {
    buildTree();
    BH = new bhdata[nNodes];
    propagateCOM(ROOT);
    checkCOM(ROOT);
    propagateBmax(ROOT);
    propagateExtent(ROOT);

    STATS(
        // space for statistics
        alloc(nleaf, nprocs);     alloc(nnode, nprocs);
        alloc(nbar, nprocs);      alloc(ndirect, nprocs);
        alloc(illl, nprocs);      alloc(illn, nprocs);
    );
    resetStats();
}

template <typename T>
void BHtree<T>::propagateCOM(int node) {

    BH[node].com.zero();
    BH[node].mass = 0;

    if( isLeaf(node) ) {
        FORALLPOINTS(i, node) {
            BH[node].com += ps[i].mass * ps[i].r;
            BH[node].mass += ps[i].mass;
        }
    }
    else {
        FORALLCHILDREN(d,node) {
            propagateCOM(d);
            BH[node].com  += BH[d].mass * BH[d].com;
            BH[node].mass += BH[d].mass;
        }
    }

    BH[node].com /= BH[node].mass;
}

template <typename T>
void BHtree<T>::checkCOM(int node) {

    vec com(0);
    T mass = 0;
    int nn = 0;
    FORALLPOINTS(i, node) {
        com += ps[i].mass * ps[i].r;
        mass += ps[i].mass;
        nn++;
    }
    com /= mass;

    const T tol =  std::numeric_limits<T>::epsilon()*10000;
    if( !CLOSE(com, BH[node].com, tol) ) {
        pprint("BHtree<T>::checkCOM error:\n");
        pprint("Node %d (%d):\n", node, isLeaf(node));
        pprint("% e  % e %4d  % e\n", com, BH[node].com, nn, com - BH[node].com);
        exit(1);
    }

    FORALLCHILDREN(d, node) checkCOM(d);
}


// Extent of a leaf is just it's bmax -- because it doesn't have daughters
template <typename T>
void BHtree<T>::propagateExtent( int node ) {
    double extent = 0;
    if ( isLeaf( node ) ) {
        extent = BH[node].Bmax;
    }
    else {
        FORALLCHILDREN ( d, node ) {
            propagateExtent( d );
            extent = std::max( extent, ( BH[d].com - BH[node].com ).norm() + BH[d].extent );
        }
    }
    extent = std::max( extent, BH[node].Bmax );
    assert( extent >= BH[node].Bmax );
    BH[node].extent = extent;
}


template <typename T>
void BHtree<T>::propagateBmax(int node) {
    T bmax = 0, b2 = 0, b3 = 0, b4 = 0;
    FORALLPOINTS(i, node) {
        T dr = (ps[i].r - BH[node].com).norm();
        bmax = std::max( bmax, dr );
        b2 += dr * dr;
        b3 += dr * dr * dr;
        b4 += dr * dr * dr * dr;
    }
    BH[node].Bmax = bmax;
    BH[node].B2 = b2;
    BH[node].B3 = b3;
    BH[node].B4 = b4;

    if( !isLeaf(node) ) FORALLCHILDREN(d,node)  propagateBmax(d);
}


//----------------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------------
// From here on, we have different traversals of the tree.
//----------------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------------

// Barnes-Hut multipole acceptance criterion
template <typename T>
inline bool BHtree<T>::BHThetaMAC( int p, T theta, int node ) {
    // (work in squares avoid the sqrt)
    double r2 = ( ps[p].r - BH[node].com ).norm2();
    double b2 = sqr(BH[node].extent);

    return ( b2 < sqr(theta) * r2 );
}

//----------------------------------------------------------------------------------------------------
// simple Barnes-Hut single-particle-at-a-time routines without assembly
//----------------------------------------------------------------------------------------------------

// apply a range of sources to a sink
template <typename T>
inline void BHtree<T>::direct_noAVX(int sink, int srcBeg, int srcEnd, vec delta) {
    for(int src=srcBeg; src<=srcEnd; src++) {
        vec dr = ps[sink].r - ps[src].r + delta;
        T ir = 1.0/sqrt(dr.norm2() + eps2);
        ps[sink].a -= ps[src].mass*ir*ir*ir*dr;
        ps[sink].pot -= ps[src].mass*ir;
    }

}

// apply one source to a sink
template <typename T>
inline void BHtree<T>::direct_noAVX(int sink, vec source, T srcmass, vec delta) {
    vec dr = ps[sink].r - source + delta;
    T ir = 1.0/sqrt(dr.norm2() + eps2);
    ps[sink].a -= srcmass*ir*ir*ir*dr;
    ps[sink].pot -= srcmass*ir;
}

template <typename T>
void BHtree<T>::bhSingle_noAVX(int sink, T theta, int node, vec delta) {

    if(isLeaf(node)) {
        // apply the sink to the node
        direct_noAVX(sink, begin(node), end(node), delta);
        STATS(illl[g] += nodeSize(node); ndirect[g] += nodeSize(node);)
    }

    else {
        // monopole BH MAC
        if( BHThetaMAC(sink, theta, node ) ) {
            // apply the node to the sink
            direct_noAVX(sink, BH[node].com, BH[node].mass, delta);
            STATS(illn[g]++;)
        }
        else
            FORALLCHILDREN(d,node) bhSingle_noAVX(sink, theta, d, delta);
    }

}

template <typename T>
void BHtree<T>::acc1_noAVX(int sink, T theta, vec delta) {
    bhSingle_noAVX(sink, theta, ROOT, delta);
    STATS(nleaf[g]++;  nbar[g]++;)
}

template <typename T>
void BHtree<T>::accAll_noAVX(T theta, vec delta) {
#pragma omp parallel for schedule(dynamic, 1)
    for(int sink=0; sink<np; sink++) {
        ps[sink].a.zero();
        ps[sink].pot = 0;
        bhSingle_noAVX(sink, theta, ROOT, delta);
    }
}

//----------------------------------------------------------------------------------------------------
// simple Barnes-Hut single-particle-at-a-time routines using AVX
//----------------------------------------------------------------------------------------------------

template <typename T>
void BHtree<T>::bhSingle(int sink, T theta, int node, vec delta) {

    if(isLeaf(node)) {
        int g = omp_get_thread_num();
        avx[g].addSources( ps, begin(node), end(node) );
    }
    else {
        if( BHThetaMAC(sink, theta, node) ) {
            int g = omp_get_thread_num();
            avx[g].addSrcPoint( BH[node].com, BH[node].mass );
        }
        else {
            FORALLCHILDREN(d,node) bhSingle(sink, theta, d, delta);

        }
    }
}

template <typename T>
void BHtree<T>::acc1(int p, T theta, vec delta) {
    int g = omp_get_thread_num();
    avx[g].resetSrcCount();
    bhSingle(p, theta, ROOT, delta);
    avx[g].compute(ps, p, p, delta, eps2);

    STATS(nleaf[g]++; illl[g]+=avx[g].nsrc; ndirect[g]+=avx[g].nsrc; nbar[g]++;)
}

template <typename T>
void BHtree<T>::accAll(T theta, vec delta) {
#pragma omp parallel for schedule(dynamic, 1)
    for(int p=0; p<np; p++) {
        ps[p].a.zero();
        ps[p].pot = 0;
        acc1(p, theta, delta);
    }
}


//----------------------------------------------------------------------------------------------------
// Recursive Subset Partition Traversal:
// Find the interaction list one a node at a time, subject to its not becoming too large.
// This allows nodes which are not leaves to be used as sinks. Obviously this will be most
// useful if the maxleafsize is not too large.
// (we should also try the variant which decays to single-particle BH when necessary)
//----------------------------------------------------------------------------------------------------

/*
  If sink is well-separated from src, add src as a multipole representation.
  else
     if src is a leaf
        add its particles to interaction list
     else
        open src
*/
template <typename T>
void BHtree<T>::bhTraverse(T sinkBmax, vec sinkcom, int srcnode, T theta, int maxmax, bool flag, int thread) {

    // node-node MAC (different from point-node MAC)
    T dr2 = (sinkcom - BH[srcnode].com).norm2();
    int WS = sqr(sinkBmax + BH[srcnode].Bmax) < sqr(theta)*dr2;

    // If srcnode is well-separated from the sinknode, add the monopole representation of srcnode
    // to the interaction list.
    // if this makes the interaction list too long, abort so that we will open the sinknode
    if( WS ) {
        if(flag) if(avx[thread].nsrc + 1 > maxmax) return;
        avx[thread].addSrcPoint(BH[srcnode].com, BH[srcnode].mass);
    }
    else {
        // Otherwise, if srcnode is a leaf, we have no recourse but to add its contents as sources.
        // Thus, we are doing point-on-point between the sink and source nodes, so this is "exact"
        // Again, if this makes the interaction list too long, go back to open sink node
        if( isLeaf(srcnode) ) {
            if(flag) if(avx[thread].nsrc + nodeSize(srcnode) > maxmax) return;
            avx[thread].addSources( ps, begin(srcnode), end(srcnode) );
        }
        // otherwise, open the srcnode
        else {
            FORALLCHILDREN(d,srcnode) bhTraverse(sinkBmax, sinkcom, d, theta, maxmax, flag, thread);
        }
    }
}


template <typename T>
void BHtree<T>::recursiveSubsetPartitionInternal(int sinknode, T theta, int maxILlength, vec delta, int level) {

    // We're going to start computing the acceleration of the contents of sinknode.
    // Clear out the interaction list (stored in avx[g])
    int g = omp_get_thread_num();
    avx[g].resetSrcCount();

    // if sinknode is a leaf, do a traverse and take whatever interaction list length it requires
    // (the "false" parameter).
    // We won't consider opening it further and doing individual BH computations for now...
    if( isLeaf(sinknode) ) {
        bhTraverse(BH[sinknode].Bmax, BH[sinknode].com, ROOT, theta, maxILlength, false, g);
        //DIRECT(g, sinknode);
        avx[g].compute(ps, begin(sinknode), end(sinknode), delta, eps2);
        STATS( nleaf[g]++; illl[g]+=avx[g].nsrc; ndirect[g]+= avx[g].nsrc*nodeSize(sinknode);
               nbar[g]+=nodeSize(sinknode); )
    }
    else {
        // Otherwise, if sinknode is not a leaf, do a traverse, but stop if the interaction list length
        // becomes too large and open the sinknode to try again
        bhTraverse(BH[sinknode].Bmax, BH[sinknode].com, ROOT, theta, maxILlength, true, g);
        if(avx[g].nsrc >= maxILlength) {
            FORALLCHILDREN(d, sinknode) {
#pragma omp task
                recursiveSubsetPartitionInternal(d, theta, maxILlength, delta, level+1);
            }
        }
        else {
            // the travese resulted in a sufficiently short interaction list length
            //DIRECT(g, sinknode);
            avx[g].compute(ps, begin(sinknode), end(sinknode), delta, eps2);
            STATS( nnode[g]++; illn[g]+=avx[g].nsrc; ndirect[g]+=avx[g].nsrc*nodeSize(sinknode);
                   nbar[g]+=nodeSize(sinknode); )
        }
    }
}

template <typename T>
void BHtree<T>::recursiveSubsetPartition(T theta, int maxILlength, vec delta) {
#pragma omp parallel
    {
#pragma omp single
        recursiveSubsetPartitionInternal(0, theta, maxILlength, delta, 1);
    }
}

//----------------------------------------------------------------------------------------------------
// Traversal which uses leaves as sink cells:
// Find the interaction list for each leaf, no matter how large.
// (we should also try the variant which decays to single-particle BH when necessary)
//----------------------------------------------------------------------------------------------------

template <typename T>
void BHtree<T>::BHsubsetsInternal(int sinknode, T theta, int maxILlength, vec delta, int level) {

    // We're going to start computing the acceleration on the contents of sinknode if it
    // is a leaf.

    // if sinknode is a leaf, do a traverse and take whatever interaction list length it requires
    // (the "false" parameter).
    // We won't consider opening it further and doing individual BH computations for now...
    if( isLeaf(sinknode) ) {
        // clear out the interaction list (stored in avx[g])
        int g = omp_get_thread_num();
        avx[g].resetSrcCount();

        bhTraverse(BH[sinknode].Bmax, BH[sinknode].com, ROOT, theta, maxILlength, false, g);
        avx[g].compute(ps, begin(sinknode), end(sinknode), delta, eps2);
        STATS( nleaf[g]++; illl[g]+=avx[g].nsrc; ndirect[g]+=avx[g].nsrc*nodeSize(sinknode);
               nbar[g]+=nodeSize(sinknode); )
    }
    else {
        // otherwide travese to next leaves
        FORALLCHILDREN(d, sinknode) {
#pragma omp task
            BHsubsetsInternal(d, theta, maxILlength, delta, level+1);
        }
    }
}

template <typename T>
void BHtree<T>::BHsubsets(T theta, int maxILlength, vec delta) {

#pragma omp parallel
    {
#pragma omp single
        BHsubsetsInternal(0, theta, maxILlength, delta, 0);
    }
}

#endif // __BHTREE_INCLUDE__
