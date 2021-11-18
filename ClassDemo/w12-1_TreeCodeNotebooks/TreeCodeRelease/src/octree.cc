/*
Time-stamp: <octree.cc on Monday, 15 November, 2021 at 23:29:41 MST (pinto)>

  basic octree in 3 dimensions

  with current version of MortonKey class, can do 42 levels
  i.e. each particle coordinate in the Morton key is represented by an unsigned 42 bit integer
  (about 12 digits of precision in scaled coordinates).

*/

#ifndef __OCTREE_INCLUDE__
#define __OCTREE_INCLUDE__

#include <numeric>
#include "inttype.h"
#include "smallvec.cc"
#include "pprint.cc"
#include "permute.cc"
#include "morton.cc"
#include "particleset.cc"

#define FORALLCHILDREN(D,node)    for(int D=tree[node].firstChild; D<tree[node].firstChild+tree[node].childCount; D++)
#define FORALLPOINTS(I,node)      for(int I=tree[node].begin; I<=tree[node].end; I++)

// for printing from Python
void printit(const char *s) {
    pprint("%s", s);
}

template <typename T>
class Octree : public Morton {
public:

private:
    using vec = SmallVec<T, 3>;
    using keyType = Morton::mortonKey;
    using ptype = Ptype<T>;

    // node of the tree:
    struct TreeNode {

        TreeNode() {
            cellCentre = {0.0,0.0,0.0};
            halfWidth = 0;
            begin = end = 0;
            parent = 0;
            firstChild = whichChild = childCount = level = 0;
        }

        vec   cellCentre;   // centre coordinate of cell
        double halfWidth;    // half-width of the cell
        uint32 begin, end;   // beginning and ending particle number (inclusive) in node
        uint32 parent;       // node of parent
        uint32 firstChild;   // node of first child
        uint16 whichChild;   // my direction from parent
        uint8  childCount;   // how many children (0 to 8)
        uint8  level;        // level in tree
    };

public:
    Octree(ParticleSet<T> &particleSet, int maxLeafSize);
    ~Octree();

    void buildTree();
    void checkTree();

    // make the morton keys & put points in Morton order
    void makeMorton();
    void passOne();
    void passTwo();
    void sortPoints();

    // ordering of particles
    void putParticlesInTreeOrder();
    void putParticlesInOriginalOrder();
    bool particlesInTreeOrder();

    // the two traversals; pass one counts the number of nodes at each level.
    void inOrderTraversalPassOne(int parentLevel, int b, int e);
    // pass two fills in the tree nodes once it has been allocated.
    void inOrderTraversalPassTwo(int level, int b, int e, int parent, vec centre, T hw);

    // checks the correctness of the tree
    void traverseCheckTree(int node, int level, vec nodecentre, T hw);

    // some geometry functions
    bool within(const vec p, vec nodecentre, T hw);
    void reportWithin(vec pos, vec nodecentre, T hw, int level);
    int internalWhereAmI(typename Morton::mortonKey mymorton, int node, int level);
    int sphereNodeIntersect(vec sinkcentre, T r2, int node);
    int sphereContainsNode(vec const centre, T const r2, int const node);

    // beginning/ending particle indices in a node
    inline int begin(const int node)  { assert( node < nNodes ); return tree[node].begin; }
    inline int end(const int node)    { assert( node < nNodes ); return tree[node].end; }

    // number of particles in node
    inline int nodeSize(const int node)  { assert(node < nNodes); return tree[node].end - tree[node].begin + 1; }

    // is particle w/ index p in node?
    bool isIn(int p, int node) { return (p>=begin(node)) && (p<=end(node)); }

    // is node a leaf?
    inline bool isLeaf(const int node)       { assert(node < nNodes); return (tree[node].childCount == 0); }

    // return coordinates of center of node
    inline vec nodeCentre(const int node)  { assert(node < nNodes); return tree[node].centre; }

    // determine the node number of the leaf cell where this position belongs
    int whereAmI(vec pos);

    static const int maxLevel = 42;
    bool verbose = false;      // print information about the tree
    bool check = false;        // perform various checks

    int np;                    // number of particles
    ParticleSet<T> &ps; // particle set

    keyType *morton;           // Morton keys
    int *pindex;               // permutation array
    int *rindex;               // reverse permutation array
    bool inTreeOrder;          // are the particles in Morton order?

    int maxLeafSize;           // parameter for tree
    int maxChildren;           // maximum number of children is 8 for three dimensions
    vec ROOTcentre;            // center of problem
    T halfWidth;               // halfwidth of problem


    TreeNode *tree;            // the tree itself
    int nLeaves, *leafSet;     // list of leaf nodes
    int nNodes, *node2leaf;    // map from leafSet index into tree node
    int avgLeafSize;
    int ROOT, ROOTlevel, nodeptr;
    int levelcount[maxLevel], levelptr[maxLevel];

    // allow printing from C++ or from Python
    void (*printn)(const char *s) = printit;

    template<typename... Args>
    void myprint(Args... args) {
        printn(spprint(args...).c_str());
    }

    void registerPrint( void (*print_callback)(const char *s) ) {
        printn = print_callback;
    }

};

template <typename T>
Octree<T>::Octree(ParticleSet<T> &particleSet, int maxLeafSize) : ps(particleSet), maxLeafSize(maxLeafSize) {
    maxChildren = 8;
    ROOT = 0;    ROOTlevel = 0;
    inTreeOrder = false;

    // pointers to internal arrays
    tree = NULL; morton = NULL; pindex = NULL; rindex = NULL; leafSet = NULL; node2leaf = NULL;
}

template <typename T>
Octree<T>::~Octree() {
    // delete all allocated memory
    if(      tree != NULL ) delete[] tree;
    if(    morton != NULL ) delete[] morton;
    if(    pindex != NULL ) delete[] pindex;
    if(    rindex != NULL ) delete[] rindex;
    if(   leafSet != NULL ) delete[] leafSet;
    if( node2leaf != NULL ) delete[] node2leaf;
}

template <typename T>
void Octree<T>::buildTree() {

    np = ps.N;

    // take care scaled particle positions must be on [0,1)^3 (i.e., must be >=0, < 1.0)
    assert( ps.hasBoundingBox() );
    ROOTcentre = ps.centre;
    halfWidth = ps.halfWidth * (1.0 + pow(2.0, -maxLevel));

    // make Morton keys
    makeMorton();

    // sort particles and keys
    sortPoints();

    // first pass to determine number of node per level in the tree
    // so that tree nodes at each level are stored consecutively
    passOne();

    // second pass to fill in tree nodes
    passTwo();

    // check tree if requested
    if(check) {
        checkTree();
        myprint("tree passes tests\n");
    }
}

// Make Morton Keys
template <typename T>
void Octree<T>::makeMorton() {

    // upper limit slightly larger so that scaled particle position are strictly < boxMax
    this->initMorton(ROOTcentre-vec(ps.halfWidth), ROOTcentre+vec(halfWidth));

    morton = new typename Morton::mortonKey[np];
    for(int i=0; i<np; i++) morton[i] = this->makeKey(ps[i].r);

}

// sort points into Morton order & create permutation arrays pindex and rindex
template <typename T>
void Octree<T>::sortPoints() {

    // permutation which takes Morton keys (hence particles) to Morton order
    pindex = new int[np];

    // fill with 0,1,2,...
    std::iota(pindex, pindex+np, 0);

    // find permutation array which sorts particles into Morton order
    // last argument is lambda function for sorting the permuation array
    // according to comparisons of the Morton keys...
    sortWithPermutation(morton, pindex, np,
                        [&](const int &a, const int &b) { return AM(morton[a], morton[b]); }
                        );

    if(check) {
        for(int i=0; i<np-1; i++) assert( AM(morton[pindex[i]], morton[pindex[i+1]]) );
    }

    // create rindex to put them back in their original order
    rindex = new int[np];
    for(int i=0; i<np; i++) rindex[pindex[i]] = i;

    // actually permute particles and their indices into Morton order
    ps.permute(pindex, np);
    // and keep the Morton keys in the same order
    inSituPermute(morton, pindex, np);
    inTreeOrder = true;
}

// run through the Morton keys once to determine space required for tree
// and number of nodes per level.
template <typename T>
void Octree<T>::passOne() {

    nLeaves = 0;
    levelcount[0] = 1;
    for(int l=1; l<maxLevel; l++) levelcount[l] = 0;

    nNodes = 1; // count the root node!
    avgLeafSize = 0;
    inOrderTraversalPassOne(ROOTlevel, 0, np-1);

    if(verbose) {
        myprint("levelcount: \n");
        for(int i=0; i<maxLevel; i++) myprint(" %d,",levelcount[i]);
        myprint("\n");

        myprint("maxLeafSize: %d\n", maxLeafSize);
        myprint("nLeaves: %d\n", nLeaves);
        myprint("<leaf size> = %e\n", (double)avgLeafSize/nLeaves);
    }

}

// run through the Morton keys a second time to fill in the tree node data.
template <typename T>
void Octree<T>::passTwo() {

    // find number of cells in tree
    int nc = 0;
    for(int l=0; l<maxLevel; l++) nc += levelcount[l];
    if(verbose) myprint("nc: %d   nNodes: %d\n", nc, nNodes);
    assert( nc == nNodes );

    // allocate space for leafSet and node2leaf index mapping
    node2leaf = new int[nNodes];
    for(int i=0;i<nNodes;i++) node2leaf[i] = -1;

    leafSet = new int[nLeaves];
    nLeaves = 0;

    // allocate space for tree
    tree = new TreeNode[nNodes]();

    // levelptr[i] points to next node at level i to be filled in
    levelptr[0] = 0;
    for(int l=1; l<maxLevel; l++) levelptr[l] = levelptr[l-1] + levelcount[l-1];

    // pointer to current node in tree
    nodeptr = ROOT;

    // set up root node
    tree[nodeptr].cellCentre = ROOTcentre;  tree[nodeptr].halfWidth = halfWidth;
    tree[nodeptr].begin = 0;                tree[nodeptr].end = np-1;
    tree[nodeptr].firstChild = 0;           tree[nodeptr].childCount = 0;
    tree[nodeptr].parent = -1;              tree[nodeptr].whichChild = 0;
    tree[nodeptr].level = 0;

    // run through the data a second time to build the tree
    inOrderTraversalPassTwo(ROOTlevel, 0, np-1, ROOT, ROOTcentre, 0.5*halfWidth);
    assert( nodeptr+1 == nNodes );
}

template <typename T>
bool Octree<T>::particlesInTreeOrder() {
    return inTreeOrder;
}

template <typename T>
void Octree<T>::putParticlesInTreeOrder() {
    if( !inTreeOrder ) {
        //inSituPermute(ps, pindex, np);
        ps.permute(pindex, np);
        inTreeOrder = true;
    }
}

template <typename T>
void Octree<T>::putParticlesInOriginalOrder() {
    if( inTreeOrder ) {
        //inSituPermute(ps, rindex, np);
        ps.permute(rindex, np);
        inTreeOrder = false;
    }
}

// first pass to count cells and leaves
template <typename T>
void Octree<T>::inOrderTraversalPassOne(int parentLevel, int b, int e) {

    // get the direction to first child
    int direction = this->levelKey(morton[b],parentLevel);

    while( (b<=e) ) {
        assert (direction<=maxChildren);

        // count number of particles in this octant
        int count = 0;
        while( b <= e ) {
            if( this->levelKey(morton[b],parentLevel) == direction ) { b++; count++; }
            else break;
        }
        assert(count>0);

        // increment number of nodes at this parentLevel
        levelcount[parentLevel+1]++;
        nNodes++;

        // if leaf node, go no deeper
        if( count <= maxLeafSize ) {
            nLeaves++;
            avgLeafSize += count;
        }
        else
            inOrderTraversalPassOne(parentLevel+1, b-count, b-1);

        // get the drection to next child
        if(b<=e) direction = this->levelKey(morton[b],parentLevel);
    }
}

// second pass to fill in tree data; same traversal as in pass one.
template <typename T>
void Octree<T>::inOrderTraversalPassTwo(int parentLevel, int b, int e, int parent, vec centre, T hw) {

    assert( parentLevel < maxLevel ); assert( e>b );
    // shouldn't be here if parent has any children
    assert( tree[parent].firstChild == 0 ); assert( tree[parent].childCount == 0 );

    // get the direction to first child
    int direction = this->levelKey(morton[b],parentLevel);

    while( b<=e ) {
        assert( direction <= maxChildren );

        // count number of particles in this octant
        int count = 0;
        while( b <= e ) {
            if( this->levelKey(morton[b],parentLevel) == direction ) { b++; count++; }
            else break;
        }
        assert( count > 0 );

        // get child node number in tree
        int child = levelptr[parentLevel+1];
        levelptr[parentLevel+1]++;

        if( tree[parent].firstChild == 0 ) {
            // first child
            assert( tree[parent].childCount == 0 );
            tree[parent].firstChild=child;
            tree[parent].childCount = 1;
        }
        else {
            assert( tree[parent].childCount <= maxChildren );
            tree[parent].childCount++;
        }

        tree[child].level = parentLevel+1;
        tree[child].parent = parent;
        tree[child].whichChild = direction;
        tree[child].begin = b-count;
        tree[child].end = b-1;
        tree[child].halfWidth = hw;
        tree[child].cellCentre = centre + hw*this->mask[direction];
        tree[child].firstChild = 0;
        tree[child].childCount = 0;

        nodeptr++;
        assert( nodeptr < nNodes );

        if( count <= maxLeafSize ) {
            // node has <= maxLeafSize particles -- make it a leaf
            leafSet[nLeaves] = child;
            node2leaf[child] = nLeaves;
            nLeaves++;
        }
        else {
            // node has too many particles -- recurse to divide into children
            inOrderTraversalPassTwo(tree[child].level, b-count, b-1, child, tree[child].cellCentre, 0.5*hw);
        }

        // more points left, get next octant to process
        if(b<=e) direction = this->levelKey(morton[b],parentLevel);
    }
}

//------------------------------------------------------------------------------------------------------------------------
// checking code
//------------------------------------------------------------------------------------------------------------------------

// is a particle within the node
template <typename T>
bool Octree<T>::within(const vec pos, vec nodecentre, T hw) {
    for(int d=0; d<3; d++)
        if( !(nodecentre[d] - hw <= pos[d] && pos[d] < nodecentre[d] + hw) ) return false;
    return true;
}

// format complaint
template <typename T>
void Octree<T>::reportWithin(vec pos, vec nodecentre, T hw, int level) {
    myprint("\n\nParticle out of place within tree:\n");

    std::string name[3] = {"x", "y", "z"};
    myprint("nodecentre: % .15e  hw: % .15e\n", nodecentre, hw);
    for(int d=0; d<3; d++)
        myprint("%s:  % .15e <= % .15e <= % .15e\n",
               name[d], nodecentre[d]-hw, pos[d], nodecentre[d]+hw);
    vec up = nodecentre+vec(hw);
    vec dn = nodecentre-vec(hw);
    this->printKey(this->makeKey(dn));
    this->printKey(this->makeKey(pos));
    this->printKey(this->makeKey(up));
    for(int i=1; i<27+level-1; i++) myprint(" "); myprint("^\n");
    myprint("level: %d\n", level);
    vec pos_scaled = pos - this->boxMin;
    up -= this->boxMin;
    dn -= this->boxMin;
    for(int d=0; d<3; d++) {
        pos_scaled[d] *= this->scale[d];
        up[d] *= this->scale[d];
        dn[d] *= this->scale[d];
        myprint("%d <= %d <= %d\n", dn[d], pos_scaled[d], up[d]);
        myprint("\n");
    }
}

// traverse the tree checking that each point is within its purported node and
// the levels, widths, and centres are correct
template <typename T>
void Octree<T>::traverseCheckTree(int node, int level, vec nodecentre, T hw) {

    assert( tree[node].level == level );
    assert( tree[node].halfWidth == hw );

    FORALLPOINTS(p, node) {
        if( !within(ps[p].r, nodecentre, hw) ) reportWithin(ps[p].r, nodecentre, hw, level);
        assert( within( ps[p].r, nodecentre, hw) );
    }

    FORALLCHILDREN(child, node) {
        vec dcentre = nodecentre;
        dcentre += 0.5*hw*Morton::mask[tree[child].whichChild];
        vec dc = tree[child].cellCentre;
        assert(dc == dcentre );
        traverseCheckTree(child, level+1, dcentre, 0.5*hw);
    }
}

// public interface for checking tree
template <typename T>
void Octree<T>::checkTree() {
    traverseCheckTree(ROOT, ROOTlevel, ROOTcentre, halfWidth);
}

//------------------------------------------------------------------------------------------------------------------------
// utility geometry functions
//------------------------------------------------------------------------------------------------------------------------

// given a MortonKey, determine the node number of the leaf cell where this key belongs
template <typename T>
int Octree<T>::internalWhereAmI(typename Morton::mortonKey mymorton, int node, int level) {

    // if a leaf, terminate with my node
    if(isLeaf(node)) return node;

    // else recurse into the correct direction
    int direction = this->key(mymorton,level);
    int child = -1;
    FORALLCHILDREN(d, node) {
        if( tree[d].wc == direction ) {
            child = internalWhereAmI(mymorton, d, level+1);
            break;
        }
    }
    if( child == -1 ) {
        myprint("internalWhereAmI: leaf cell doesn't exist in tree\n");
        assert( child >= 0 );
    }

    return child;
}

// return the index of the node containing pos,  giving an error if no leaf node contains this pos
template <typename T>
int Octree<T>::whereAmI(vec pos) {
    typename Morton::mortonKey mymorton = this->makeKey(pos);
    return internalWhereAmI(mymorton, ROOT, ROOTlevel);
}

// sphere-cube intersection; returns true if any part of node is within the sphere
template <typename T>
int Octree<T>::sphereNodeIntersect(vec centre, T r2, int node) {
    assert(node < nNodes);
    vec nodecentre = tree[node].centre;

    // fast rejection: sphere-sphere intersection
    T r = sqrt(r2);
    if( (centre-nodecentre).norm2() > sqr( sqrt(3)*tree[node].halfWidth + r ) ) return 0;

    // they might intersect -- do the more expensive exact test
    T hw = tree[node].halfWidth;
    vec nmin = nodecentre - vec(hw);
    vec nmax = nodecentre + vec(hw);

    T mindist2 = 0;
    for(int d=0; d<3; d++) {
        if(      centre[d] < nmin[d] ) mindist2 += sqr( centre[d] - nmin[d] );
        else if( centre[d] > nmax[d] ) mindist2 += sqr( centre[d] - nmax[d] );
    }
    return mindist2 <= r2;
}

template <typename T>
int Octree<T>::sphereContainsNode(vec const centre, T const r2, int const node) {
    assert(node < nNodes);
    vec nodecentre = tree[node].centre;

    // fast rejection: sphere-sphere intersection
    T r = sqrt(r2);
    if( (centre-nodecentre).norm2() > sqr( sqrt(3)*tree[node].halfWidth + r ) ) return 0;

    // find vertex furthest from centre
    vec vtx;
    for(int d=0; d<3; d++) vtx[d] = nodecentre[d] + ( nodecentre[d] > centre[d] ? 1 : -1 ) * tree[node].halfWidth;
    return (vtx-centre).norm2() <= r2;
}

#endif // __OCTREE_INCLUDE__


#ifdef TEST_OCTREE

#include "STimer.cc"

int main() {

    using Real = double;
    using vec = SmallVec<Real, 3>;

    srand48(93847593);

    int np = 1000000;
    ParticleSet<Real> pset(np);

    for(int i=0; i<np; i++) {
        vec r = {drand48(), drand48(), drand48()};
        r = 5 * r - vec(2.5, 5.0, 3.0);
        vec v = {drand48(), drand48(), drand48()};
        vec a(0.0);
        double pot = 0;
        double m = 1/np;

        // Why?
        //pset[i] = {r, v, a, pot, m, i};
        PtypeRef<Real> Q =  {r, v, a, pot, m, i};
        pset[i] = Q;
    }
    pset.computeBoundingBox();

    STimer treeTimer;

    int maxLeafSize = 32;
    Octree<Real> O = Octree<double>(pset, maxLeafSize);
    O.verbose = true;
    O.check = true;


    treeTimer.CStart();
    O.buildTree();
    treeTimer.Stop();
    treeTimer.Report("time to build tree:");
}
#endif // TEST_OCTREE
