import torch

def flux(rho, gamma):
    return gamma * rho * (1. - rho)

def d_flux(rho, gamma):
    return gamma*(1. - 2.*rho)

def Rusanov_Flux(rho, gamma):
    left = rho[:-1]
    right = rho[1:]

    s = torch.tensor([max(abs(d_flux(rho[j], gamma)), abs(d_flux(rho[j+1], gamma))) for j in range(len(rho)-1)]) 

    return 0.5*(flux(left, gamma) + flux(right, gamma)) - 0.5 * s * (right - left)