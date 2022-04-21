# used pretty much everywhere
import numpy as np
# used to enable function overloading
from multimethod import multimeta



class SoA(object, metaclass=multimeta):
    _pref  = "particle_"
    _pos   = "position_"
    _vel   = "vel"

    _id    = "id"
    _cpu   = "cpu"
    _id_0  = "id_0"
    _cpu_0 = "cpu_0"
    _id_1  = "id_1"
    _cpu_1 = "cpu_1"


    def __init__(self, data, copy_v=False, copy_id=False):

        # Flag which fields are copied
        self.contains_v  = copy_v
        self.contains_id = copy_id

        # Copy positions
        str_pos = self._pref + self._pos
        self.px = np.array(data["all", str_pos + "x"])
        self.py = np.array(data["all", str_pos + "y"])
        self.pz = np.array(data["all", str_pos + "z"])

        # Copy velocities
        if copy_v:
            str_vel = self._pref + self._vel
            self.vx = np.array(data["all", str_vel + "x"])
            self.vy = np.array(data["all", str_vel + "y"])
            self.vz = np.array(data["all", str_vel + "z"])

        # Copy particle ID
        if copy_id:
            self.id    = np.array(data["all", self._pref + self._id],    dtype=int)
            self.cpu   = np.array(data["all", self._pref + self._cpu],   dtype=int)
            self.id_0  = np.array(data["all", self._pref + self._id_0],  dtype=int)
            self.cpu_0 = np.array(data["all", self._pref + self._cpu_0], dtype=int)
            self.id_1  = np.array(data["all", self._pref + self._id_1],  dtype=int)
            self.cpu_1 = np.array(data["all", self._pref + self._cpu_1], dtype=int)



    def __print_pos(self):
        return str(self.px) + "," + str(self.py) + "," + str(self.pz)


    def __print_vel(self):
        return str(self.vx) + "," + str(self.vy) + "," + str(self.vz)


    def __print_id(self):
        return "id:"   + str(self.id)   + ", cpu:"   + str(self.cpu) + ", " + \
               "id_0:" + str(self.id_0) + ", cpu_0:" + str(self.cpu_0) + ", " + \
               "id_1:" + str(self.id_1) + ", cpu_1:" + str(self.cpu_1)



    def __str__(self):
        str_rep = "{pos:"  + self.__print_pos()
        if self.contains_v:
            str_rep += "; vel:" + self.__print_vel()
        if self.contains_id:
            str_rep += "; " + self.__print_id()
        str_rep += "}"

        return str_rep


    def __repr__(self):
        return str(self)


    def pos(self, index: int):
        return self.px[index], self.py[index], self.pz[index]

    def pos(self, start: int, stop: int):
        return self.px[start:stop], self.py[start:stop], self.pz[start:stop]

    def pos(self):
        return self.px, self.py, self.pz


    def vel(self, index: int):
        if self.contains_v:
            return self.vx[index], self.vy[index], self.vz[index]
        else:
            return None

    def vel(self, start:int , stop: int):
        if self.contains_v:
            return self.vx[start:stop], self.vy[start:stop], self.vz[start:stop]
        else:
            return None

    def vel(self):
        if self.contains_v:
            return self.vx, self.vy, self.vz
        else:
            return None


    def pid(self, index: int):
        if self.contains_id:
            return self.id[index],   self.cpu[index], \
                   self.id_0[index], self.cpu_0[index], \
                   self.id_1[index], self.cpu_1[index]
        else:
            return None

    def pid(self, start: int, stop: int):
        if self.contains_id:
            return self.id[start:stop],   self.cpu[start:stop], \
                   self.id_0[start:stop], self.cpu_0[start:stop], \
                   self.id_1[start:stop], self.cpu_1[start:stop]
        else:
            return None


    def pid(self):
        if self.contains_id:
            return self.id, self.cpu, self.id_0, self.cpu_0, self.id_1, self.cpu_1
        else:
            return None



class Particle(object):
    def __init__(self, px, py, pz, **kwargs):
        self.pos = np.array([px, py, pz])

        self.contains_vel = False
        if "vel" in kwargs.keys():
            self.vel = np.array(kwargs["vel"])
            self.contains_vel = True

        self.contains_id = False
        if "id" in kwargs.keys():
            self.id    = kwargs["id"][0]
            self.cpu   = kwargs["id"][1]
            self.id_0  = kwargs["id"][2]
            self.cpu_0 = kwargs["id"][3]
            self.id_1  = kwargs["id"][4]
            self.cpu_1 = kwargs["id"][5]
            self.contains_id = True


    def __str__(self):
        str_rep = "P(" + str(self.pos)
        if self.contains_vel:
            str_rep += ", " + str(self.vel)
        if self.contains_id:
            str_rep += ", " + str(self.id) +   ", " + str(self.cpu)   + \
                       ", " + str(self.id_0) + ", " + str(self.cpu_0) + \
                       ", " + str(self.id_1) + ", " + str(self.cpu_1)
        str_rep += ")"

        return str_rep


    def __repr__(self):
        return str(self)



class AoS(object):
    def __init__(self, amrex_data, copy_v=False, copy_id=False):
        self.particles = list()
        soa = SoA(amrex_data, copy_v=copy_v, copy_id=copy_id);

        for i, elt in enumerate(zip(*soa.pos())):
            if soa.contains_v and soa.contains_id:
                self.particles.append(Particle(* elt, vel=soa.vel(i), id=soa.pid(i)))
            elif soa.contains_v:
                self.particles.append(Particle(* elt, vel=soa.vel(i)))
            elif soa.contains_id:
                self.particles.append(Particle(* elt, id=soa.pid(i)))
            else:
                self.particles.append(Particle(* elt))
