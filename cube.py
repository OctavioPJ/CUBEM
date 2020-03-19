# -*- coding: utf-8 -*-
"""
Created on Wed Feb  6 11:33:21 2019

@author: CNEA
"""
import sys
if 'D:\\' not in sys.path: sys.path.append('D:\\')
from PyCAR.PyCIT.lector_mest_posta import LoadCDP, get_reactivity_pro as grp, extraer_bet, AddableDict
from PyCAR.PyCIT.FT import MeshFlux
import matplotlib.pyplot as plt
import numpy as np
import os
import time

BLACKABSORBER = 1999
DEFAULT_STATE = 1
DEFAULT_BURNUP = 0


class CUBE_MAP(object):
    _wdir = 'C:\\CUBEM\\'

    def MapMaterial(self, meshmaterial, **kwargs):
        if meshmaterial == BLACKABSORBER:
            return self._NULL
        return self._FISS[1][0.0]


class Nu_Fission_Map(CUBE_MAP):
    def __init__(self, **kwargs):
        if 'wdir' in kwargs:
            setattr(self, '_wdir', kwargs['wdir'])

        wdir = getattr(self, '_wdir')
        self._FISS = LoadCDP(wdir + 'ec18cs1gf.cdp')
        self._NULL = {'XS': np.zeros(self._FISS[1][0.0]['XS'].shape), \
                      'SM': np.zeros(self._FISS[1][0.0]['SM'].shape)}


class Kinetic_Map(CUBE_MAP):
    def __init__(self, **kwargs):
        if 'wdir' in kwargs:
            setattr(self, '_wdir', kwargs['wdir'])
        wdir = getattr(self, '_wdir')

        self._FISS = extraer_bet(wdir + 'ec18cs1gf.cdo')

        if 'NPRC' in kwargs and kwargs['NPRC'] == 1:
            print('CONDENSED')
            condensed_kin = {'Beta': np.array([self._FISS[DEFAULT_STATE][DEFAULT_BURNUP]['Beta'][-1]])}
            condensed_kin.update({'Decays':
                                      np.array([condensed_kin['Beta'] / (
                                              self._FISS[DEFAULT_STATE][DEFAULT_BURNUP]['Beta'][:-1] /
                                              self._FISS[DEFAULT_STATE][DEFAULT_BURNUP]['Decays']).sum()])})
            condensed_kin.update({'Neutron Velocity': self._FISS[DEFAULT_STATE][DEFAULT_BURNUP]['Neutron Velocity']})
            self._FISS[DEFAULT_STATE][DEFAULT_BURNUP].update(condensed_kin)

        self._NULL = {KinParam: np.zeros(self._FISS[DEFAULT_STATE][DEFAULT_BURNUP][KinParam].shape) \
                      for KinParam in ['Beta', 'Decays', 'Neutron Velocity']}

    def get_NPRC(self, *args, **kwargs):
        if args:
            return self._FISS[DEFAULT_STATE][DEFAULT_BURNUP][args[0]].shape[0]
        return self._FISS[DEFAULT_STATE][DEFAULT_BURNUP]['Decays'].shape[0]


if __name__ == '__main__':

    from citvappru.SourceCAREM2_1_1 import Geometry, PrecCalc
    import re

    os.chdir('C:\\CUBEM\\')
    FILENAME = 'cube.cii'
    ROOTFILE = FILENAME[:-4]
    DATABASE = ROOTFILE + '_eq.cdb'

    Geom = Geometry(FILENAME.replace('.cii', '0.cii'))

    NuFis = Nu_Fission_Map()
    KM = Kinetic_Map(NPRC=1)

    _NOG = 1

    NPRC = KM.get_NPRC()

    Flux = MeshFlux(DATABASE + '.meshflux', *Geom.Cantidad_de_Nodos(), _NOG)
    Nx, Ny, Nz = Flux.shape[1:4]

    BetaM = np.empty((1, Nx, Ny, Nz, NPRC))
    LambdaM = np.empty((1, Nx, Ny, Nz, NPRC))
    NuFisM = np.empty((1, Nx, Ny, Nz, _NOG))
    VelocityM = np.empty((1, Nx, Ny, Nz, _NOG))

    state = 0

    Vmesh = Geom.Vmesh()

    #    react = grp(ROOTFILE+'.plt')
    react = [173.0]
    keff = 1/(1-react[-1]/100/1000)
    # keff = 1.0015

    for _x in range(Nx):
        for _y in range(Ny):
            for _z in range(Nz):
                meshmaterial = Geom.sc5[_x][_y][_z]

                KinParam = KM.MapMaterial(meshmaterial, NPRC=NPRC)

                BetaM[state][_x][_y][_z][:] = KinParam['Beta']

                LambdaM[state][_x][_y][_z][:] = KinParam['Decays']

                VelocityM[state][_x][_y][_z][:] = KinParam['Neutron Velocity']

                NuFisM[state][_x][_y][_z][:] = NuFis.MapMaterial(meshmaterial)['XS'][:, 3]

    TimeSteps = np.logspace(2, 5, num=10, dtype=int)[3:]

    for N in TimeSteps:

        C0 = {state: {}}
        for nx in range(Flux.shape[1]):
            C0[state][nx] = {}
            for ny in range(Flux.shape[2]):
                C0[state][nx][ny] = {}
                for nz in range(Flux.shape[3]):
                    FluxL = [Flux[state, nx, ny, nz, group] for group in range(Flux.shape[-1])]
                    NuFisL = [NuFisM[state][nx][ny][nz][group] / keff for group in range(Flux.shape[-1])]
                    Nu_FluxM = [NuFisL[group] * FluxL[group] * Vmesh for group in range(Flux.shape[-1])]

                    Bet_k = BetaM[state][nx][ny][nz]
                    Lamb_k = LambdaM[state][nx][ny][nz]
                    Nu_Flux = sum(Nu_FluxM)

                    C0[state][nx][ny][nz] = [Bet_k[prec] * Nu_Flux / Lamb_k[prec]
                                                 if Lamb_k[prec] != 0 else 0.0
                                             for prec in range(NPRC)]
        alfi = AddableDict(C0).as_array()

        p = re.compile(r'POWER\(WATTS\)\s+([0-9]\.[0-9]{5}E\+[0-9]{2})')
        Q = {}
        EQUILIBRIUM = True
        chi_g = [1.0]

        tfinal = 3.0

        DATABASE = ROOTFILE + 'S.cdb'
        Times = []
        powers = []
        precursors = []
        source = []
        dt = tfinal/N
        FAILED_COMPLETE_TEST = False
        DERIVATIVE = True
        for n_step in range(N+1):
            t = n_step * dt
            if EQUILIBRIUM:
                C_t_1 = C0.copy()
                C_t = C0.copy()
                Flux_t_1 = Flux
            else:
                TFlux = MeshFlux(DATABASE + '.meshflux', Nx, Ny, Nz, _NOG)
                C_t = PrecCalc(C_t_1, TFlux, NuFisM, Vmesh, dt
                               , LambdaM=LambdaM
                               , BetaM=BetaM)

                C_t_1 = C_t.copy()
                Flux_t_1 = TFlux

            for group in range(Flux.shape[-1]):
                Q[group] = {}
                for state in range(Flux.shape[0]):
                    Q[group][state] = {}
                    for nx in range(Flux.shape[1]):
                        Q[group][state][nx] = {}
                        for ny in range(Flux.shape[2]):
                            Q[group][state][nx][ny] = {}
                            for nz in range(Flux.shape[3]):
                                _Lmk = LambdaM[state][nx][ny][nz]
                                _C = C_t[state][nx][ny][nz]

                                _invV = VelocityM[state, nx, ny, nz, group]
                                T_1Flux = Flux_t_1[state, nx, ny, nz, group]

                                Q[group][state][nx][ny][nz] = \
                                    chi_g[group] * sum([_Lmk[prc] * _C[prc] for prc in range(NPRC)]) \
                                    + _invV / dt * T_1Flux * Vmesh * DERIVATIVE

            with open('source.dat', 'w') as fod:
                for group in range(_NOG):
                    for state in range(Flux.shape[0]):
                        for nz in range(Flux.shape[3]):
                            for ny in range(Flux.shape[2]):
                                for nx in range(Flux.shape[1]):
                                    fod.write('{:15.7E}'.format(Q[group][state][nx][ny][nz]))
                    fod.write('\n')

            OsArgs = (FILENAME.replace('.cii', 'S.cii'),
                      _NOG, Nx, Ny, Nz, EQUILIBRIUM,
                      *['{:12.5E}'.format(_d) for _d in [dt * DERIVATIVE, t]])

            os.system('ExSource.bat ' + ' '.join(map(str, OsArgs)))
            # time.sleep(0.1)
            fid = open(ROOTFILE + 'S.cio', 'r')
            powers.append(float(next(p.finditer(fid.read())).groups()[0]))
            print(powers[-1])
            precursors.append(AddableDict(C_t).as_array().sum())
            source.append(AddableDict(Q).as_array().sum())
            fid.close()
            if EQUILIBRIUM:
                EQUILIBRIUM = False
            #                FAILED_COMPLETE_TEST = True
        #                break
            Times.append(t)
            assert len(powers) == len(Times)
        assert t == 3.0, 'Tiempo {}'.format(t)

        # if FAILED_COMPLETE_TEST:
            # Times = np.array([k * dt for k in range(len(powers))])

        Normalized_Pow = np.array(powers) / 100E+06
        NParray = Normalized_Pow.reshape(len(Normalized_Pow), 1)
        NTimes = np.array(Times).reshape((len(Times), 1))
        np.savetxt('Evolution{}.txt'.format(dt), np.r_['1', NTimes, NParray])
        # plt.plot(Times, Normalized_Pow, 'o', Times, precursors / precursors[0], 'o', Times, source / source[0], 'o')
