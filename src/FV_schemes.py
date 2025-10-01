import torch


def flux(rho, gamma):
    return gamma * rho * (1. - rho)

def d_flux(rho, gamma):
    return gamma*(1. - 2.*rho)


def Rusanov_Flux(rho, gamma):
    left = rho[:-1]
    right = rho[1:]

    a = torch.tensor([max(abs(d_flux(rho[j], gamma)), abs(d_flux(rho[j+1], gamma))) for j in range(len(rho)-1)]) 

    return 0.5*(flux(left, gamma) + flux(right, gamma)) - 0.5 * a * (right - left)

def Rusanov_Flux_2(left, right, gamma):
    a = torch.max(torch.abs(d_flux(left, gamma)), torch.abs(d_flux(right, gamma)))
    
    return 0.5*(flux(left, gamma) + flux(right, gamma)) - 0.5 * a * (right - left)





# WICH OF THE FOLLOWING FUNCTION ARE NEEDED?? PROBABLY NONE
sigma = torch.tensor(0.5)

def fmax(gamma):
    return flux(sigma, gamma)

def D(rho, gamma):
    if rho <= sigma:
        return flux(rho, gamma)
    else:
        return fmax(gamma)

def S(rho, gamma):
    if rho <= sigma:
        return fmax(gamma)
    else:
        return flux(rho, gamma)
