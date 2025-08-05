import numpy as np
import matplotlib.pyplot as plt

class LightCurveMultiple(object):

    def __init__(self, ticksinper = 100, depth = 0, duration = 0, noise = .001, numper = 1):
        """
        Light Curve with transit details for one period

        Args:
        ticksinper (int): number of timesteps in one period
        length (integer): number of total timesteps
        depth (0 < float < 1.0): depth of simulated transit
        duration (0 < float < 1.0): time length of simulated transit
        noise (float): normalized noise of flux
        location (0 < float < 1.0): central location of transit
        flux (array of floats): normalized flux of light curve
        numper (int): number of periods

        """
        self.ticksinper = ticksinper
        self.length = self.ticksinper * numper
        if depth == 0:
            self.depth = np.random.uniform(.01,.20)
        else:
            self.depth = np.random.normal(loc = depth, scale = .03, size = 1)
        if duration == 0:
            self.duration = np.random.randint(.02 * ticksinper, .15 * ticksinper)
        else:
            self.duration = duration * ticksinper
        self.noise = noise
        self.location = np.random.randint(0,ticksinper)
        self.flux = np.ones(self.length) + np.random.normal(loc = 0, scale = self.noise, size = self.length)
        self.per = 0
        self.numper = numper
    
    def new_transit(self):
        self.depth = np.random.normal(loc=self.depth, scale = 0.001, size = 1)
        self.location = self.location + self.ticksinper
        
    def plot_transit(self, phase_flag = False):
        """ 
        Subtracts transit from the flux and plots the resulting lightcurve
        """
        lowerbound = int(self.location - self.duration/2)
        upperbound = int(self.location + self.duration/2)
        #print("FIRST")
        #print(lowerbound,upperbound,self.location)
        if lowerbound < 0:
            self.ub = self.length + lowerbound
            self.lb = upperbound
            self.flux[0:self.lb] -= self.depth
            self.flux[self.ub:self.length] -= self.depth
            overflow_flag = True
        elif upperbound > self.length:
            self.lb = int(upperbound - self.length)
            self.ub = int(lowerbound)
            self.flux[0:self.lb] -= self.depth
            self.flux[self.ub:self.length] -= self.depth
            self.overflow_flag = True
        else:
            self.ub = int(upperbound)
            self.lb = int(lowerbound) 
            self.flux[self.lb:self.ub] -= self.depth
            self.overflow_flag = False
        self.plot(phase_flag=True)

    def plot(self, phase_flag = False):
        
        if phase_flag == True:
            if self.overflow_flag == True:
                if self.per == 0:
                    plt.plot(np.arange(self.ticksinper*self.per+self.lb+1, self.ticksinper * (1 + self.per))/self.ticksinper % 1, self.flux[self.ticksinper*self.per+self.lb+1:self.ticksinper * (1 + self.per)], label = "Lightcurve", color = "deepskyblue")
                    plt.plot(np.arange(0,self.lb+1)/self.ticksinper,self.flux[0:self.lb+1] % 1, label = "Transit", color = "navy", zorder = 3)
                else:
                    plt.plot(np.arange(self.ticksinper*self.per+self.lb+1, self.ticksinper * (1 + self.per))/self.ticksinper % 1, self.flux[self.ticksinper*self.per+self.lb+1:self.ticksinper * (1 + self.per)], label = "Lightcurve", color = "deepskyblue")
                    plt.plot(np.arange(0,self.lb+1)/self.ticksinper,self.flux[0:self.lb+1] % 1, color = "navy", zorder = 3)
                plt.plot(np.arange(self.ub-1,self.length)/self.ticksinper, self.flux[self.ub-1:self.length] % 1, color = "navy", zorder = 3)
                plt.xlabel("Period")
                plt.ylabel("Normalized Flux")
                plt.title("Generated Light Curve")
                plt.legend()            
            else:
                if self.per == 0:
                    plt.plot(np.arange(self.ticksinper*self.per, self.ticksinper * (1 + self.per))/self.ticksinper % 1, self.flux[self.ticksinper*self.per:self.ticksinper * (1 + self.per)], label = "Lightcurve", color = "deepskyblue")
                    plt.plot(np.arange(self.lb-1,self.ub+1)/self.ticksinper % 1,self.flux[self.lb-1:self.ub+1], label = "Transit", color = "navy", zorder = 3)
                else:
                    plt.plot(np.arange(self.ticksinper*self.per, self.ticksinper * (1 + self.per))/self.ticksinper % 1, self.flux[self.ticksinper*self.per:self.ticksinper * (1 + self.per)], color = "deepskyblue")
                    plt.plot(np.arange(self.lb-1,self.ub+1)/self.ticksinper % 1,self.flux[self.lb-1:self.ub+1], color = "navy", zorder = 3)
                plt.xlabel("Period")
                plt.ylabel("Normalized Flux")
                plt.title("Generated Light Curve")
                plt.legend()
        else:

            if self.overflow_flag == True:
                if self.per == 0:
                    plt.plot(np.arange(self.ticksinper*self.per+self.lb+1, self.ticksinper * (1 + self.per))/self.ticksinper, self.flux[self.ticksinper*self.per+self.lb+1:self.ticksinper * (1 + self.per)], label = "Lightcurve", color = "deepskyblue")
                    plt.plot(np.arange(0,self.lb+1)/self.ticksinper,self.flux[0:self.lb+1], label = "Transit", color = "navy", zorder = 3)
                else:
                    plt.plot(np.arange(self.ticksinper*self.per+self.lb+1, self.ticksinper * (1 + self.per))/self.ticksinper, self.flux[self.ticksinper*self.per+self.lb+1:self.ticksinper * (1 + self.per)], color = "deepskyblue")
                    plt.plot(np.arange(0,self.lb+1)/self.ticksinper,self.flux[0:self.lb+1], color = "navy", zorder = 3)                    
                plt.plot(np.arange(self.ub-1,self.length)/self.ticksinper, self.flux[self.ub-1:self.length], color = "navy", zorder = 3)
                plt.xlabel("Period")
                plt.ylabel("Normalized Flux")
                plt.title("Generated Light Curve")
                plt.legend()            
            else:
                if self.per == 0:
                    plt.plot(np.arange(self.ticksinper*self.per, self.ticksinper * (1 + self.per))/self.ticksinper, self.flux[self.ticksinper*self.per:self.ticksinper * (1 + self.per)], label = "Lightcurve", color = "deepskyblue")
                    plt.plot(np.arange(self.lb-1,self.ub+1)/self.ticksinper,self.flux[self.lb-1:self.ub+1], label = "Transit", color = "navy", zorder = 3)
                else:
                    plt.plot(np.arange(self.ticksinper*self.per, self.ticksinper * (1 + self.per))/self.ticksinper, self.flux[self.ticksinper*self.per:self.ticksinper * (1 + self.per)], color = "deepskyblue")
                    plt.plot(np.arange(self.lb-1,self.ub+1)/self.ticksinper,self.flux[self.lb-1:self.ub+1], color = "navy", zorder = 3)
                plt.xlabel("Period")
                plt.ylabel("Normalized Flux")
                plt.title("Generated Light Curve")
                plt.legend()
        self.per += 1

