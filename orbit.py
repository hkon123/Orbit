import math
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import numpy as np
global G
G = 6.67408*10**-11

class Moon(object):

    def __init__(self, mass, radius):
        self.mass = mass
        self.radius = radius
        self.position = np.array([self.radius,0])
        self.velocity = None
        self.accs = None
        self.revaccs = np.array([(-G)*(self.mass/self.radius**2),0])


    def Mvelocity(self,i):
         return math.sqrt(self.velocity[i,0]**2+self.velocity[i,1]**2)

class Planet(object):

    def __init__(self, mass):
        self.mass = mass
        self.moons = []
        self.position = np.array([0,0])
        self.velocity = np.array([0,0])
        self.energies = []

    def Pvelocity(self,i):
        return math.sqrt(self.velocity[i,0]**2+self.velocity[i,1]**2)

    def addMoon(self,moons):
        for i in moons:
            self.moons.append(i)
            i.accs = np.array([(-G)*(self.mass/i.radius**2),0])
            i.velocity = np.array([0,math.sqrt(((G)*self.mass)/i.radius)])

    def distance(self,i):
        diff = self.position[i]-self.moons[0].position[i]
        length = math.sqrt(diff[0]**2+diff[1]**2)
        return length
    def direction(self,i):
        direct = self.position[i]-self.moons[0].position[i]
        return direct

    def energy(self):
        initial = 0.5*self.moons[0].mass*self.moons[0].Mvelocity(0)**2
        quarter = 0.5*self.moons[0].mass*self.moons[0].Mvelocity(self.steps/4)**2 + 0.5*self.mass*self.Pvelocity(self.steps/4)**2
        half = 0.5*self.moons[0].mass*self.moons[0].Mvelocity(self.steps/2)**2 + 0.5*self.mass*self.Pvelocity(self.steps/2)**2
        threequarter = 0.5*self.moons[0].mass*self.moons[0].Mvelocity(3*self.steps/4)**2 + 0.5*self.mass*self.Pvelocity(3*self.steps/4)**2
        final = 0.5*self.moons[0].mass*self.moons[0].Mvelocity(self.steps-1)**2 + 0.5*self.mass*self.Pvelocity(self.steps-1)**2

        self.energies = [initial,quarter,half,threequarter,final]
        

    def sim(self,stepLength,steps):
        self.steps = steps
        for k in range(0,steps):
            for i in self.moons:
                if k ==0:
                    i.velocity = np.vstack((i.velocity,i.velocity + i.accs*stepLength))
                    i.position = np.vstack((i.position,i.position + i.velocity[k+1]*stepLength))
                    
                    self.velocity = np.vstack((self.velocity, self.velocity + i.revaccs*stepLength))
                    self.position = np.vstack((self.position, self.position + self.velocity[k+1]*stepLength))
                    
                    i.accs = (self.direction(k+1)/self.distance(k+1))*(G)*(self.mass/i.radius**2)
                    i.revaccs = (self.direction(k+1)/self.distance(k+1))*(-G)*(i.mass/i.radius**2)
                else:
                    i.velocity = np.vstack((i.velocity,i.velocity[k] + i.accs*stepLength))
                    i.position = np.vstack((i.position,i.position[k] + i.velocity[k+1]*stepLength))
                    
                    self.velocity = np.vstack((self.velocity, self.velocity[k] + i.revaccs*stepLength))
                    self.position = np.vstack((self.position, self.position[k] + self.velocity[k+1]*stepLength))
                    
                    i.accs = (self.direction(k+1)/self.distance(k+1))*(G)*(self.mass/i.radius**2)
                    i.revaccs = (self.direction(k+1)/self.distance(k+1))*(-G)*(i.mass/i.radius**2)
                                                                                                      
    def init(self):
        # initialiser for animator
        return self.patches

    def animate(self, i):
        self.patches[0].center = (self.moons[0].position[i,0], self.moons[0].position[i,1])
        self.patches[1].center =(self.position[i,0],self.position[i,1])
        return self.patches

    def run(self):
        # create plot elements
        fig = plt.figure()
        ax = plt.axes()

        # create circle of radius 0.1 centred at initial position and add to axes
        self.patches = []
        
        self.patches.append(plt.Circle((self.moons[0].position[0,0], self.moons[0].position[0,1]), 22333, color = 'g', animated = True))
        self.patches.append(plt.Circle((0,0),6792000, color = 'r', animated = True))
        for i in range(0,len(self.patches)):
            ax.add_patch(self.patches[i])

        # set up the axes
        ax.axis('scaled')
        ax.set_xlim(-(self.moons[0].radius+10**6), self.moons[0].radius+10**6)
        ax.set_ylim(-(self.moons[0].radius+10**6), self.moons[0].radius+10**6)
        #ax.set_xlabel('x (rads)')
        #ax.set_ylabel('sin(x)')

        # create the animator
        anim = FuncAnimation(fig, self.animate, init_func = self.init, frames = self.steps, repeat = False, interval = 1, blit = True)

        # show the plot
        plt.show()


def main():
    phobos = Moon(1.06*10**16,9.3773*10**6)
    mars = Planet(6.4185*10**23)
    mars.addMoon([phobos])
    mars.sim(50,3000)
    mars.run()
    #diff = []
    #for i in range(0,len(phobos.position)):
        #a = math.sqrt(phobos.position[i,0]**2+phobos.position[i,1]**2)
        #diff.append(a-9.3773*10**6)
    #print diff
    mars.energy()
    print(mars.energies)


main()
