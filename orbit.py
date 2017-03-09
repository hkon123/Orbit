import math
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import numpy as np
global G
G = 6.67408*10**-11

class Moon(object):

    def __init__(self, mass, radius,size, color):
        self.color = color
        self.size = size
        self.mass = mass
        self.radius = radius
        self.position = np.array([self.radius,0])
        self.velocity = None
        self.accs = None
        self.revaccs = np.array([(-G)*(self.mass/self.radius**2),0])
        self.period = None


    def Mvelocity(self,i):
         return math.sqrt(self.velocity[i,0]**2+self.velocity[i,1]**2)

    def mdistance(self,moon1,moon2,i):
        diff = moon1.position[i]-moons2.position[i]
        length = math.sqrt(diff[0]**2+diff[1]**2)
        return length
    def mdirection(self,i):
        direct = moon1.position[i]-moons2.position[i]
        return direct

class Planet(object):

    def __init__(self, mass, size,color):
        self.color = color
        self.size = size
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

    def distance(self,moon,i):
        diff = self.position[i]-moon.position[i]
        length = math.sqrt(diff[0]**2+diff[1]**2)
        return length
    def direction(self, moon, i):
        direct = self.position[i]-moon.position[i]
        return direct

    def energy(self):
        initial =0
        quarter =0
        half = 0
        threequarter = 0
        final = 0
        for i in range(0,len(self.moons)):
            initial += 0.5*self.moons[i].mass*self.moons[i].Mvelocity(0)**2
        for i in range(0,len(self.moons)):
            quarter += 0.5*self.moons[i].mass*self.moons[i].Mvelocity(self.steps/4)**2
        quarter += 0.5*self.mass*self.Pvelocity(self.steps/4)**2
        for i in range(0,len(self.moons)):
            half += 0.5*self.moons[i].mass*self.moons[i].Mvelocity(self.steps/2)**2
        half += 0.5*self.mass*self.Pvelocity(self.steps/2)**2
        for i in range(0, len(self.moons)):
            threequarter += 0.5*self.moons[i].mass*self.moons[i].Mvelocity(3*self.steps/4)**2
        threequarter += 0.5*self.mass*self.Pvelocity(3*self.steps/4)**2
        for i in range(0,len(self.moons)):
            final += 0.5*self.moons[i].mass*self.moons[i].Mvelocity(self.steps-1)**2
        final += 0.5*self.mass*self.Pvelocity(self.steps-1)**2

        self.energies = [initial,quarter,half,threequarter,final]
        

    def sim(self,stepLength,steps):
        self.steps = steps
        self.stepLength = stepLength
        for k in range(0,steps):
            for i in self.moons:
                if k ==0:
                    i.velocity = np.vstack((i.velocity,i.velocity + i.accs*stepLength))
                    i.position = np.vstack((i.position,i.position + i.velocity[k+1]*stepLength))
                    
                    self.velocity = np.vstack((self.velocity, self.velocity + i.revaccs*stepLength))
                    self.position = np.vstack((self.position, self.position + self.velocity[k+1]*stepLength))

                    if len(self.moons)<1:
                        for j in self.moons:
                            if j != i:
                                accs = (mdirection(i,j,k)/mdistance(i,j,k))*(-G)*(j.mass/mdistance(i,j,k)**2)
                                i.velocity = np.vstack((i.velocity,i.velocity + accs*stepLength))
                                i.position = np.vstack((i.position,i.position + i.velocity[k+1]*stepLength))
                    
                    #i.accs = (self.direction(k+1)/self.distance(k+1))*(G)*(self.mass/self.distance(k+1)**2)
                    #i.revaccs = (self.direction(k+1)/self.distance(k+1))*(-G)*(i.mass/self.distance(k+1)**2)
                else:

                    i.accs = (self.direction(i,k)/self.distance(i,k))*(G)*(self.mass/self.distance(i,k)**2)
                    i.revaccs = (self.direction(i,k)/self.distance(i,k))*(-G)*(i.mass/self.distance(i,k)**2)
                    
                    i.velocity = np.vstack((i.velocity,i.velocity[k] + i.accs*stepLength))
                    i.position = np.vstack((i.position,i.position[k] + i.velocity[k+1]*stepLength))
                    
                    self.velocity = np.vstack((self.velocity, self.velocity[k] + i.revaccs*stepLength))
                    self.position = np.vstack((self.position, self.position[k] + self.velocity[k+1]*stepLength))

                    if len(self.moons)<1:
                        for j in self.moons:
                            if j != i:
                                accs = (mdirection(i,j,k)/mdistance(i,j,k))*(-G)*(j.mass/mdistance(i,j,k)**2)
                                i.velocity = np.vstack((i.velocity,i.velocity + accs*stepLength))
                                i.position = np.vstack((i.position,i.position + i.velocity[k+1]*stepLength))

                    if (i.position[k,1]==0) or(i.position[k-1,1]<0 and i.position[k,1]>0):
                        if i.period == None:
                            i.period = float(k*self.stepLength)/(60*60*24)
                    
                    
                                                                                                      
    def init(self):
        # initialiser for animator
        return self.patches

    def animate(self, i):
        for j in range(0,len(self.moons)):
            self.patches[j].center = (self.moons[j].position[i,0], self.moons[j].position[i,1])
        self.patches[len(self.moons)].center =(self.position[i,0],self.position[i,1])
        return self.patches

    def run(self):
        # create plot elements
        fig = plt.figure()
        ax = plt.axes()

        # create circle of radius 0.1 centred at initial position and add to axes
        self.patches = []
        for j in range(0,len(self.moons)):
            self.patches.append(plt.Circle((self.moons[j].position[0,0], self.moons[j].position[0,1]), self.moons[j].size, color = self.moons[j].color, animated = True))
        self.patches.append(plt.Circle((0,0),self.size, color = self.color, animated = True))
        for i in range(0,len(self.patches)):
            ax.add_patch(self.patches[i])

        # set up the axes
        ax.axis('scaled')
        ax.set_xlim(-(self.moons[0].radius+10**7), self.moons[0].radius+10**7)
        ax.set_ylim(-(self.moons[0].radius+10**7), self.moons[0].radius+10**7)
        #ax.set_xlabel('x (rads)')
        #ax.set_ylabel('sin(x)')

        # create the animator
        anim = FuncAnimation(fig, self.animate, init_func = self.init, frames = self.steps, repeat = False, interval = 1, blit = True)

        # show the plot
        plt.show()


def main():



    
    phobos = Moon(1.06*10**16,9.3773*10**6, 223330,'g')
    deimos = Moon(1.80*10**15,23.463*10**6,124000,'b')
    mars = Planet(6.4185*10**23,6792000,'r')
    #mars.addMoon([phobos])
    mars.addMoon([deimos])
    mars.addMoon([phobos])
    mars.sim(50,3000)
    mars.run()
    #diff = []
    #for i in range(0,len(phobos.position)):
    #    a = math.sqrt(phobos.position[i,0]**2+phobos.position[i,1]**2)
    #    diff.append(a-9.3773*10**6)
    #print diff
    mars.energy()
    print(deimos.period)


main()
