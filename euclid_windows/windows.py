import numpy as np
from scipy import integrate

class Windows:
    def __init__(self,
        zmax=2.5,
        nbin=3,
        nz=101,
        dz=0.001,
        sigma_ph=0.05,
        bintype="equipopulated",

    ):

        self.zmax=zmax
        self.nbin=nbin
        self.nz=nz
        self.dz=dz
        self.sigma_ph=sigma_ph
        self.bintype=bintype   

    def get_bins(self):
        n_tot, error = integrate.quad(galaxy_distribution, 0, self.zmax)

        if self.bintype == "equipopulated":
            z_bin_edge = np.zeros(self.nbin + 1, "float64")
            total_count = 0.0
            for Bin in range(self.nbin - 1):
                bin_count = 0.0
                z = z_bin_edge[Bin]
                while bin_count <= (n_tot - total_count) / (self.nbin - Bin):
                    gd_1 = galaxy_distribution(z)
                    gd_2 = galaxy_distribution(z + self.dz)
                    bin_count += 0.5 * (gd_1 + gd_2) * self.dz
                    z += self.dz
                z_bin_edge[Bin + 1] = z
                total_count += bin_count
            z_bin_edge[self.nbin] = self.zmax
        elif self.bintype == "equispaced":
            z_bin_edge = np.linspace(0, self.zmax, self.nbin + 1)

        return z_bin_edge

    def get_distributions(self):
        self.z_bin_edge = self.get_bins()
        zeta = np.linspace(0, self.zmax, num=self.nz)

        eta_z = np.zeros((self.nz, self.nbin), "float64")
        gal = galaxy_distribution(zeta)
        for Bin in range(self.nbin):
            low = self.z_bin_edge[Bin]
            hig = self.z_bin_edge[Bin + 1]
            for inz in range(self.nz):
                z = zeta[inz]
                integrand = gal * photo_z_distribution(z, zeta, self.sigma_ph)
                integrand = np.array(
                    [
                        elem if low <= zeta[index] <= hig else 0
                        for index, elem in enumerate(integrand)
                    ]
                )
                eta_z[inz, Bin] = integrate.trapz(
                    integrand,
                    zeta,
                )

        eta_norm = np.zeros(self.nbin, "float64")
        for Bin in range(self.nbin):
            eta_norm[Bin] = np.sum(
                0.5 * (eta_z[1:, Bin] + eta_z[:-1, Bin]) * (zeta[1:] - zeta[:-1])
            )

        for Bin in range(self.nbin):
            eta_z[:, Bin] /= eta_norm[Bin]

        self.zeta = zeta
        self.eta_norm = eta_norm
        self.eta_z = eta_z
        self.bias = np.sqrt(1+(self.z_bin_edge[1:]+self.z_bin_edge[:-1])/2)

        return


def galaxy_distribution(z, zmean = 0.9):
    """
    Galaxy distribution returns the function D(z) from the notes

    """
    z0 = zmean/np.sqrt(2.0)

    galaxy_dist = z**2*np.exp(-(z/z0)**(1.5))

    return galaxy_dist

def photo_z_distribution(z, zph, sigma_ph):
    """
    Photo z distribution
    """
    photo_z_dist = np.exp(-0.5*(
        (z-zph)/sigma_ph/(1.+z))**2)/sigma_ph/(1.+z)/np.sqrt(
        2.*np.pi)

    return photo_z_dist
