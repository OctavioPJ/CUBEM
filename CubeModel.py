# -*- coding: utf-8 -*-
"""
Created on Wed Feb  6 11:33:21 2019

@author: CNEA
"""
import sys

if 'D:\\' not in sys.path: sys.path.append('D:\\')
from PyCAR.PyCIT.lector_mest_posta import LoadCDP, get_reactivity_pro as grp, extraer_bet, AddableDict
from PyCAR.PyCIT.FT import MeshFlux
from CAREM.Citvap import CitvapSourceModel
from CAREM.Poinki import EvInterface
import numpy as np
import os
from citvappru.SourceCAREM2_1_1 import Geometry
import re
import subprocess
import time
from abc import ABC, abstractmethod
BLACKABSORBER = 1999
DEFAULT_STATE = 1
DEFAULT_BURNUP = 0


class CUBE_MAP(ABC):
    _wdir = 'C:\\CUBEM\\'

    @property
    @abstractmethod
    def _NULL(self):
        return

    @property
    @abstractmethod
    def _FISS(self):
        return

    def MapMaterial(self, meshmaterial, **kwargs):
        if meshmaterial == BLACKABSORBER:
            return self._NULL
        return self._FISS[1][0.0]


class Nu_Fission_Map(CUBE_MAP):
    @property
    def _NULL(self):
        return self._NULL

    @property
    def _FISS(self):
        return self._FISS

    @_NULL.setter
    def _NULL(self, NewDict):
        assert isinstance(NewDict, dict)
        self._NULL = NewDict
        return

    @_FISS.setter
    def _FISS(self, NewValue):
        assert isinstance(NewValue, dict)
        self._FISS = NewValue

    def __init__(self, **kwargs):
        if 'wdir' in kwargs:
            setattr(self, '_wdir', kwargs['wdir'])

        wdir = getattr(self, '_wdir')
        self._FISS = LoadCDP(wdir + 'ec18cs1gf.cdp')
        self._NULL = {'XS': np.zeros(self._FISS[1][0.0]['XS'].shape),
                      'SM': np.zeros(self._FISS[1][0.0]['SM'].shape)}


class Kinetic_Map(CUBE_MAP):
    @property
    def _NULL(self):
        return self._NULL

    @property
    def _FISS(self):
        return self._FISS

    @_NULL.setter
    def _NULL(self, NewDict):
        assert isinstance(NewDict, dict)
        self._NULL = NewDict
        return

    @_FISS.setter
    def _FISS(self, NewValue):
        assert isinstance(NewValue, dict)
        self._FISS = NewValue

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

        self._NULL = {KinParam: np.zeros(self._FISS[DEFAULT_STATE][DEFAULT_BURNUP][KinParam].shape)
                      for KinParam in ['Beta', 'Decays', 'Neutron Velocity']}

    def get_NPRC(self, *args, **kwargs):
        if args:
            return self._FISS[DEFAULT_STATE][DEFAULT_BURNUP][args[0]].shape[0]
        return self._FISS[DEFAULT_STATE][DEFAULT_BURNUP]['Decays'].shape[0]


class CubeModel(EvInterface):
    @property
    def Ct_1(self):
        return AddableDict(self._Ct_1).as_array()

    @property
    def C0(self):
        return AddableDict(self._C0).as_array()

    @property
    def Ct(self):
        return AddableDict(self._Ct).as_array()

    @property
    def Q(self):
        return AddableDict(self._Q).as_array()

    def __init__(self, state=0, nprc=1,
                 Equilibrium=True, file='cubeS.cii',
                 UseDerivative=False):

        self.t = 0.0
        self.dt = 0.0
        self.state = state
        self.chi_g = [1.0]
        self.__file = file
        self.Equilibrium = Equilibrium
        self.nprc = nprc

        self.__rootfile = self.__file.replace('.cii', '')
        self.__database = self.__file.replace('.cii', '.cdb')
        self.__p_parser = re.compile(r'POWER\(WATTS\) +([0-9]\.[0-9]{5}E\+[0-9]{2})')
        self.Pt = 100
        self._Q = {}

        self._C0 = self._Ct = self._Ct_1 = {self.state: {}}

        self.__UseDerivative = UseDerivative
        self.__datfile = r'source.dat'

        self._init_citvap_model()
        self._init_mesh_parameters()
        self._init_parameters()
        self._init_flux()
        return

    def _init_citvap_model(self):
        self._CV = CitvapSourceModel(
            file=self.__file,
            sec4='cube.sc4',
            sec5='cube.sc5',
            mat='cube.mat',
            sec26='cube.sc26',
            GeometricType='XYZ')

        self._CV.PrintMeshMap()
        self._CV.Calculate()
        return

    def _init_flux(self):
        self.Flux_t = self.Flux_t_1 = \
            MeshFlux(self.__database.replace('S.cdb', '_eq.cdb.meshflux')
                                               , self.Nx, self.Ny, self.Nz, self._NOG)
        return

    def _init_parameters(self):
        self.NuFis = Nu_Fission_Map()
        self.KM = Kinetic_Map(NPRC=self.nprc)
        self._NOG = 1
        self.BetaM = np.empty((1, self.Nx, self.Ny, self.Nz, self.nprc))
        self.LambdaM = np.empty((1, self.Nx, self.Ny, self.Nz, self.nprc))
        self.NuFisM = np.empty((1, self.Nx, self.Ny, self.Nz, self._NOG))
        self.VelocityM = np.empty((1, self.Nx, self.Ny, self.Nz, self._NOG))

        for _x in range(self.Nx):
            for _y in range(self.Ny):
                for _z in range(self.Nz):
                    meshmaterial = self.Geom.sc5[_x][_y][_z]
                    KinParam = self.KM.MapMaterial(meshmaterial, NPRC=self.nprc)
                    self.BetaM[self.state][_x][_y][_z][:] = KinParam['Beta']
                    self.LambdaM[self.state][_x][_y][_z][:] = KinParam['Decays']
                    self.VelocityM[self.state][_x][_y][_z][:] = KinParam['Neutron Velocity']
                    self.NuFisM[self.state][_x][_y][_z][:] = self.NuFis.MapMaterial(meshmaterial)['XS'][:, 3]

        self.react = self.__R0 = 173.0E-5
        self.keff = 1 / (1 - self.__R0)
        return 0

    def _init_mesh_parameters(self):
        self.Geom = Geometry(self.__file.replace('S.cii', '0.cii'))  # HARDCODEADO el *S.cii
        self.Vmesh = self.Geom.Vmesh()
        self.Nx, self.Ny, self.Nz = self.Geom.Cantidad_de_Nodos()
        return

    def equilibrium_precrs(self, *args, **kwargs):
        if not self._C0[self.state]:
            for nx in range(self.Flux_t.shape[1]):
                self._C0[self.state][nx] = {}
                for ny in range(self.Flux_t.shape[2]):
                    self._C0[self.state][nx][ny] = {}
                    for nz in range(self.Flux_t.shape[3]):
                        FluxL = [self.Flux_t[self.state, nx, ny, nz, group]
                                 for group in range(self.Flux_t.shape[-1])]
                        NuFisL = [self.NuFisM[self.state][nx][ny][nz][group] / self.keff
                                  for group in range(self.Flux_t.shape[-1])]
                        Nu_FluxM = [NuFisL[group] * FluxL[group] * self.Geom.Vmesh()
                                    for group in range(self.Flux_t.shape[-1])]
                        Bet_k = self.BetaM[self.state][nx][ny][nz]
                        Lamb_k = self.LambdaM[self.state][nx][ny][nz]

                        self._C0[self.state][nx][ny][nz] = [Bet_k[prec] / Lamb_k[prec] * sum(Nu_FluxM)
                                                            if Lamb_k[prec] != 0 else 0.0 for prec in range(self.nprc)]
        return

    def calculate_precrs(self, dt):
        for nx in range(self.Nx):
            # self._Ct[self.state][nx] = {}
            for ny in range(self.Ny):
                # self._Ct[self.state][nx][ny] = {}
                for nz in range(self.Nz):
                    FluxL = [self.Flux_t[self.state, nx, ny, nz, group] for group in range(self._NOG)]
                    NuFisL = [self.NuFisM[self.state][nx][ny][nz][group] for group in range(self._NOG)]

                    Bet_k = self.BetaM[self.state][nx][ny][nz]
                    Lamb_k = self.LambdaM[self.state][nx][ny][nz]
                    _c_t_1 = self._Ct_1[self.state][nx][ny][nz]

                    Nu_Flux = [NuFisL[group] * FluxL[group] * self.Geom.Vmesh()
                               for group in range(self._NOG)]

                    self._Ct[self.state][nx][ny][nz] = \
                        [_c_t_1[prec] * np.exp(-Lamb_k[prec] * dt) +
                         Bet_k[prec] / Lamb_k[prec] * sum(Nu_Flux) * (1 - np.exp(-Lamb_k[prec] * dt))
                         if Lamb_k[prec] != 0 else 0.0 for prec in range(self.nprc)]
        return

    def evolve(self, dt, *args, **kwargs):
        if self.Equilibrium:
            self.equilibrium_precrs()
            self._Ct = self._C0.copy()
            self._Ct_1 = self._C0.copy()
        else:
            self._Ct_1 = self._Ct.copy()
            self.calculate_precrs(dt)
            self.t += dt

        if self.__UseDerivative:
            self.dt = dt

        self.calculate_source(dt)
        self.source_to_file()
        self.__send_to_lu_17__()
        self.__xsu_mod__()
        self.run()
        self.__get_power__()

        if self.Equilibrium:
            self.Equilibrium = False

        self.Flux_t_1 = self.Flux_t.copy()
        self.Flux_t = MeshFlux(self.__database+'.meshflux', self.Nx, self.Ny, self.Nz, self._NOG)
        return

    def source_to_file(self):
        with open(self.__datfile, 'w') as fod:
            for group in range(self._NOG):
                for nz in range(self.Nz):
                    for ny in range(self.Ny):
                        for nx in range(self.Nx):
                            fod.write('{:15.7E}'.format(self._Q[group][self.state][nx][ny][nz]))
                fod.write('\n')
        return

    def __send_to_lu_17__(self):
        subprocess.run(['copy', '/Y', self.__datfile, 'Source_To_LU17'], shell=True)
        os.chdir('C:\\CUBEM\\Source_To_LU17')
        msg = subprocess.run(['main.exe', *map(str, [self._NOG, self.Nx, self.Ny, self.Nz])]
                             , shell=True, capture_output=True)
        msg2 = subprocess.run(['copy', '/Y', 'fort.17', '..\\' + self.__rootfile + '.src'],
                              shell=True, capture_output=True)
        os.chdir('C:\\CUBEM')
        return msg

    def __xsu_mod__(self):
        XSU_FILE = '{}'.format(self.__rootfile + '.XSU')
        try:
            self._CV.run(executable='pre_cit.exe')
        except AssertionError:
            time.sleep(2.0)
            self._CV.run(executable='pre_cit.exe')

        subprocess.run(['move', XSU_FILE, 'XSU_MOD\\'], shell=True)
        os.chdir('C:\\CUBEM\\XSU_MOD')
        if self.__UseDerivative:
            time_step = self.dt
        else:
            time_step = 0.0
        msg = subprocess.run(['neutyp', '{}'.format(XSU_FILE), str(self.Equilibrium),
                              *map(lambda a: '{:12.5E}'.format(round(a, ndigits=5)), [time_step, self.t])],
                             shell=True, capture_output=True)
        assert msg.stderr in ['', b''], 'Error en la sobreescritura del archivo ROOT.XSU'
        os.chdir('C:\\CUBEM')
        msg2 = subprocess.run(['move', '/Y', 'C:\\CUBEM\\XSU_MOD\\{}'.format(XSU_FILE), XSU_FILE],
                              shell=True, capture_output=True)
        assert msg2.stderr in ['', b''], 'Error en el movimiento del archivo ROOT.XSU'
        return msg

    def __get_power__(self, file=None):
        if not file:
            fid = open(self.__rootfile + '.cio', 'r')
        else:
            fid = open(file, 'r')
        self.Pt = float(self.__p_parser.findall(fid.read())[0])
        fid.close()
        return

    def run(self):
        try:
            self._CV.run(executable='cittrs.exe')
        except AssertionError:
            time.sleep(2.0)
            self._CV.run(executable='cittrs.exe')
        subprocess.run(['caremdb', '-opt:export', '-val:meshflux', self.__database])
        return

    def calculate_source(self, dt):
        for group in range(self.Flux_t.shape[-1]):
            self._Q[group] = {self.state: {}}
            for nx in range(self.Flux_t.shape[1]):
                self._Q[group][self.state][nx] = {}
                for ny in range(self.Flux_t.shape[2]):
                    self._Q[group][self.state][nx][ny] = {}
                    for nz in range(self.Flux_t.shape[3]):
                        _Lmk = self.LambdaM[self.state][nx][ny][nz]
                        _C = self._Ct[self.state][nx][ny][nz]

                        inv_V = self.VelocityM[self.state, nx, ny, nz, group]

                        self._Q[group][self.state][nx][ny][nz] = \
                            self.chi_g[group] * sum([_Lmk[prc] * _C[prc] for prc in range(self.nprc)])
                        if self.__UseDerivative:
                            self._Q[group][self.state][nx][ny][nz] +=\
                                inv_V / dt * self.Flux_t_1[self.state, nx, ny, nz, group] * self.Vmesh
        return

    pass  # CubeModel


if __name__ == '__main__':
    tfinal = 3.0
    N = 100
    UseDerivative = True
    # for N in np.logspace(1, 5, dtype=int):
    for N in np.logspace(3, 5, dtype=int, num=10)[2:]:  # N = 4641
        Cube = CubeModel(UseDerivative=UseDerivative)
        dt = tfinal/N
        Times = []
        Powers = []
        for n in range(N+1):
            Cube.evolve(dt)
            print(Cube.t, Cube.Pt)
            Powers.append(Cube.Pt)
            Times.append(Cube.t)

        # assert Cube.t == 3.0
        np.savetxt('evol_{}_{}.txt'.format(N, UseDerivative),
                   np.r_['0', np.array([Times]), np.array([Powers])].T)
