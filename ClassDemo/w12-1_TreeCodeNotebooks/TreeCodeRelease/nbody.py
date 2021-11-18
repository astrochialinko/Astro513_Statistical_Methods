import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d
import numpy as np
from treeclasses import *
from timerclass import Timer

class Hamiltonian:
    def __init__(self, theta: float, epsSmooth:float, maxLeafSize:int, maxSrcLen:int):
        self.theta = theta
        self.epsSmooth = epsSmooth
        self.maxLeafSize = maxLeafSize
        self.maxSrcLen = maxSrcLen

        self.verbose = False
        self.check = False
        self.getStats = False

    def getAcceleration(self, ps: ParticleSet):

        ps.computeBoundingBox()
        center, halfWidth = ps.getBoundingBox()

        BH = BHTree(ps, self.maxLeafSize, self.epsSmooth, self.maxSrcLen)

        BH.makeTree(self.verbose, self.check, self.getStats)

        BH.BHsubsets(theta, ps.N)
        #BH.accAll(theta)

        BH.putParticlesInOriginalOrder()  # so we know which particles are which galaxy

        BH.free()                         # important to free memory on the C side...

    def positionEquation(self, ps: ParticleSet, newvel: bool) -> np.ndarray :
        return ps.v

    def momentumEquation(self, ps: ParticleSet, newacc: bool) -> np.ndarray :
        if newacc: self.getAcceleration(ps)
        return ps.a


class State:
    def __init__(self, time: float, ps: ParticleSet, hamilton: Hamiltonian ):
        self.time = time
        self.step = 0
        self.ps = ps

    def kick(self, h, hamilton, recalculate):
        # Update velocity
        self.ps.v += h * hamilton.momentumEquation(self.ps, recalculate)
        return self

    def drift(self, h, hamilton, recalculate):
        # Update positions
        self.ps.r += h * hamilton.positionEquation(self.ps, recalculate)
        return self

def KDK(dt: float, hamilton: Hamiltonian, s: State):
    s = s.kick(dt/2, hamilton, False).drift(dt, hamilton, False).kick(dt/2, hamilton, True)
    s.time += dt
    s.step += 1
    return s

def getEnergy(ps: ParticleSet):
    K = 0.5*np.sum(ps.mass * np.sum(ps.v*ps.v, axis=1))
    P = 0.5*np.sum(ps.mass * ps.pot)
    mom = np.sum(ps.mass[:,np.newaxis] * ps.v)

    return K, P, mom

with open("dubinski.tab") as f:
    data = np.loadtxt(f)

# create a particle set
PS = ParticleSet()

# allocate storage
N = data.shape[0]
PS.reserve(N)

# copy in the data
PS.mass[:] = data[:,0]
PS.r[:,:]  = data[:,1:4]
PS.v[:,:]  = data[:,4:]

# the data consists of sets of particles for the disks, bulges, and halos in each galaxy
Ndisk = 16384
Nbulge = 8192
Nhalo = 16384
Ngal = Ndisk+Nbulge+Nhalo
assert( Ngal*2 == N )

disk1 = np.s_[:Ndisk+Nbulge]
disk2 = np.s_[Ngal:Ngal+Ndisk+Nbulge]

PS.computeBoundingBox()
centre, halfWidth = PS.getBoundingBox()

fig = plt.figure()
ax = plt.axes(projection='3d')
c, hw = PS.getBoundingBox()
ax.set_xlim(c[0]-hw,c[0]+hw)
ax.set_ylim(c[0]-hw,c[0]+hw)
ax.set_zlim(c[0]-hw,c[0]+hw)
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('z')

plot1, = ax.plot3D(PS.r[disk1,0], PS.r[disk1,1], PS.r[disk1,2], 'b.', alpha=0.05, markersize=0.5)
plot2, = ax.plot3D(PS.r[disk2,0], PS.r[disk2,1], PS.r[disk2,2], 'r.', alpha=0.05, markersize=0.5)
steptxt = ax.text2D(0.2,1.0, f"step:    0   time: {0.0:8.2e}\nerror: E: {0:.2e}  P: {0:.2e}", transform=ax.transAxes)

plt.show(block=False)

def updatePlot(S, timer):
    ps = S.ps
    plot1.set_data(ps.r[disk1,0], ps.r[disk1,1]);  plot1.set_3d_properties(ps.r[disk1,2])
    plot2.set_data(ps.r[disk2,0], ps.r[disk2,1]);  plot2.set_3d_properties(ps.r[disk2,2])

    K, P, mom = getEnergy(PS)
    DErel = (K+P-E0)/E0
    DP = np.max(np.abs(mom-mom0))
    line = f"step: {S.step:4d}   time: {S.time:8.2e}\nerror: E: {DErel:.2e}  P: {DP:.2e}\nrate: {ps.N/timer.elapsed():.2e}"
    steptxt.set_text(line)

    # this updates the animation w/o stealing the window focus
    plt.gcf().canvas.draw_idle(); plt.gcf().canvas.start_event_loop(0.001)

maxLeafSize = 16
maxSrcLen = PS.N
theta = 0.75
epsSmooth = 0.025
dt = 0.1
skip = 2

H = Hamiltonian(theta, epsSmooth, maxLeafSize, maxSrcLen)
S = State(0.0, PS, H)

S = KDK(0.0, H, S) # take a step w/ dt=0 just to get initial potential
K0, P0, mom0 = getEnergy(PS)
E0 = K0+P0


def runit():
    global S
    timer = Timer()

    while(S.time<1.0):

        timer.cstart()
        S = KDK(dt, H, S) # take a KDK step
        timer.stop()

        if S.step%skip==0:
            updatePlot(S, timer)

    updatePlot(S, timer)


# profile code
# run snakeviz on nbody.prof output
if __name__ == "__main__":
    import cProfile
    import pstats

    with cProfile.Profile() as pr:
        runit()

    stats = pstats.Stats(pr)
    stats.sort_stats('tottime')
    stats.dump_stats("nbody.pprof")
