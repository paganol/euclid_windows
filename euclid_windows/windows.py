import numpy as np
from scipy import integrate
from typing import Union


class Windows:
    def __init__(
        self,
        zmin=0.001,
        zmax=2.501,
        nbin=10,
        dz=0.001,
        cb=1.0,
        zb=0,
        sigmab=0.05,
        c0=1.0,
        z0=0.1,
        sigma0=0.05,
        fout=0.1,
        bintype: Union[str, np.ndarray] = np.array([0.001, 0.418, 0.56, 0.678, 0.789, 0.9, 1.019, 1.155, 1.324, 1.576, 2.5]),
    ):

        self.dz = dz

        self.cb = cb
        self.zb = zb
        self.sigmab = sigmab
        self.c0 = c0
        self.z0 = z0
        self.sigma0 = sigma0
        self.fout = fout

        self.bintype = bintype

        if type(bintype) == np.ndarray: 
            self.nbin = len(bintype)-1
            self.zmin = bintype[0]
            self.zmax = bintype[-1]+self.dz # this adjust the last element of self.zeta
        else:
            self.nbin = nbin
            self.zmin = zmin
            self.zmax = zmax

        self.zeta = np.arange(self.zmin, self.zmax, self.dz)
        self.nz = len(self.zeta)


    def get_bins(self):
        n_tot, error = integrate.quad(galaxy_distribution, self.zmin, self.zmax)

        if any(self.bintype) == "equipopulated":
            z_bin_edge = np.zeros(self.nbin + 1, "float64") + self.zmin
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
        elif any(self.bintype) == "equispaced":
            z_bin_edge = np.linspace(self.zmin, self.zmax, self.nbin + 1)
        else:
            z_bin_edge = self.bintype

        return z_bin_edge


    def get_distributions(self):
        self.z_bin_edge = self.get_bins()

        eta_z = np.zeros((self.nz, self.nbin), "float64")
        gal = galaxy_distribution(self.zeta)
        for Bin in range(self.nbin):
            low = self.z_bin_edge[Bin]
            hig = self.z_bin_edge[Bin + 1]
            for inz in range(self.nz):
                z = self.zeta[inz]
                integrand = gal * photo_z_distribution(
                    z,
                    self.zeta,
                    self.cb,
                    self.zb,
                    self.sigmab,
                    self.c0,
                    self.z0,
                    self.sigma0,
                    self.fout,
                )
                integrand = np.array(
                    [
                        elem if low <= self.zeta[index] <= hig else 0
                        for index, elem in enumerate(integrand)
                    ]
                )
                eta_z[inz, Bin] = integrate.trapz(
                    integrand,
                    self.zeta,
                )

#        eta_norm = np.zeros(self.nbin, "float64")
#        for Bin in range(self.nbin):
#            eta_norm[Bin] = np.sum(
#                0.5 * (eta_z[1:, Bin] + eta_z[:-1, Bin]) * (self.zeta[1:] - self.zeta[:-1])
#            )
#
#        for Bin in range(self.nbin):
#            eta_z[:, Bin] /= eta_norm[Bin]

        self.eta_z = eta_z
        self.bias = np.sqrt(1 + (self.z_bin_edge[1:] + self.z_bin_edge[:-1]) / 2)

        return


def galaxy_distribution(z, zmean=0.9):
    """
    Galaxy distribution returns the function D(z) from the notes

    """
    z0 = zmean / np.sqrt(2.0)

    galaxy_dist = (z / z0) ** 2 * np.exp(-((z / z0) ** (1.5)))

    return galaxy_dist


# def photo_z_distribution(z, zph, sigma_ph):
#    """
#    Photo z distribution
#    """
#    photo_z_dist = np.exp(-0.5*(
#        (z-zph)/sigma_ph/(1.+z))**2)/sigma_ph/(1.+z)/np.sqrt(
#        2.*np.pi)
#
#    return photo_z_dist


def photo_z_distribution(
    z, zph, cb=1.0, zb=0, sb=0.05, c0=1.0, z0=0.1, s0=0.05, fout=0.1
):
    """
    Photo z distribution
    Eq. 115 and Tab. 5 of 1910.09273
    """
    sq2pi = np.sqrt(2 * np.pi)

    arg_exp_zb = -0.5 * (z - cb * zph - zb) ** 2 / (sb * (1 + z)) ** 2
    arg_exp_z0 = -0.5 * (z - c0 * zph - zb) ** 2 / (s0 * (1 + z)) ** 2

    photo_z_dist = (1 - fout) / sq2pi / sb / (1 + z) * np.exp(
        arg_exp_zb
    ) + fout / sq2pi / s0 / (1 + z) * np.exp(arg_exp_z0)

    return photo_z_dist
