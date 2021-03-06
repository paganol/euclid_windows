import numpy as np
from scipy import integrate
from typing import Union
from camb.sources import SplinedSourceWindow


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
        bintype: Union[str, np.ndarray] = np.array(
            [0.001, 0.418, 0.56, 0.678, 0.789, 0.9, 1.019, 1.155, 1.324, 1.576, 2.5]
        ),
        normalize=True,
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
            self.nbin = len(bintype) - 1
            self.zmin = bintype[0]
            self.zmax = (
                bintype[-1] + self.dz
            )  # this adjust the last element of self.zeta
        else:
            self.nbin = nbin
            self.zmin = zmin
            self.zmax = zmax

        self.zeta = np.arange(self.zmin, self.zmax, self.dz)
        self.nz = len(self.zeta)

        self.normalize = normalize

    def get_bins(self):
        self.gal_tot, error = integrate.quad(galaxy_distribution, self.zmin, self.zmax)

        if self.bintype == "equipopulated":
            z_bin_edge = np.zeros(self.nbin + 1, "float64") + self.zmin
            total_count = 0.0
            for ibin in range(self.nbin - 1):
                bin_count = 0.0
                z = z_bin_edge[ibin]
                while bin_count <= (self.gal_tot - total_count) / (self.nbin - ibin):
                    gd_1 = galaxy_distribution(z)
                    gd_2 = galaxy_distribution(z + self.dz)
                    bin_count += 0.5 * (gd_1 + gd_2) * self.dz
                    z += self.dz
                z_bin_edge[ibin + 1] = z
                total_count += bin_count
            z_bin_edge[self.nbin] = self.zmax
        elif self.bintype == "equispaced":
            z_bin_edge = np.linspace(self.zmin, self.zmax, self.nbin + 1)
        else:
            z_bin_edge = self.bintype

        self.z_bin_edge = z_bin_edge

        return

    def get_distributions(self):
        self.get_bins()

        eta_z = np.zeros((self.nbin, self.nz), "float64")
        self.gal_dist = galaxy_distribution(self.zeta)
        
        phz_dist = photo_z_distribution(
            np.array([self.zeta,]* self.nz).T,
            np.array([self.zeta,]* self.nz),
            self.cb,
            self.zb,
            self.sigmab,
            self.c0,
            self.z0,
            self.sigma0,
            self.fout,
        )

        for ibin in range(self.nbin):
            low = self.z_bin_edge[ibin]
            hig = self.z_bin_edge[ibin + 1]
            weight = np.zeros_like(self.zeta)
            weight[np.where((self.zeta >= low) & (self.zeta <= hig))] = 1.0

            eta_z[ibin, :] = self.gal_dist * integrate.trapz(
                phz_dist * weight,
                self.zeta,
                axis=1,
            )

        if self.normalize:
            eta_norm = np.zeros(self.nbin, "float64")
            for ibin in range(self.nbin):
                eta_norm[ibin] = integrate.trapz(
                    eta_z[ibin, :],
                    self.zeta,
                )

            for ibin in range(self.nbin):
                eta_z[ibin, :] /= eta_norm[ibin]

            self.eta_norm = eta_norm

        self.eta_z = eta_z
        bias = np.sqrt(1 + (self.z_bin_edge[1:] + self.z_bin_edge[:-1]) / 2)

        self.bias = np.zeros_like(self.zeta)
        for ibin in range(self.nbin):
            self.bias[
                np.where(
                    (self.zeta >= self.z_bin_edge[ibin])
                    & (self.zeta <= self.z_bin_edge[ibin + 1])
                )
            ] = bias[ibin]

        return

    def get_camb_distributions(self):
        sources = []
        try:
            self.eta_z
        except:
            print("Run before get_distributions")
            exit()

        for ibin in range(self.nbin):
            sources.append(
                SplinedSourceWindow(bias_z=self.bias, z=self.zeta, W=self.eta_z[ibin])
            )

        return sources


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

    return (1 - fout) / np.sqrt(2 * np.pi) / sb / (1 + z) * np.exp(
        -0.5 * (z - cb * zph - zb) ** 2 / (sb * (1 + z)) ** 2
    ) + fout / np.sqrt(2 * np.pi) / s0 / (1 + z) * np.exp(
        -0.5 * (z - c0 * zph - z0) ** 2 / (s0 * (1 + z)) ** 2
    )
