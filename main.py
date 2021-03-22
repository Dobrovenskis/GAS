import random
import numpy as np
import math as m

N = 100
epsilon = 1
sigma = 1
mass = 1
dt = 0.001  #sqrt(m/e) * sigma
N_steps = 10000
t = 0
n = 0.001
x_cube = (N/n)**(1/3) #sigma
V_cube = x_cube**3
V_sum_ball = V_cube/2
speed_max = 0


class atom:
    def __init__(self, pos, speed):
        self.pos = pos
        self.speed = speed
        self.pos_before = np.array([0, 0, 0], dtype = float)
        self.pos_before = pos - speed*dt
        self.a = np.array([0, 0, 0], dtype = float)

    def count_a(self, list_pos):
        F_answ = np.array([0, 0, 0], dtype = float)
        a = np.array([0, 0, 0], dtype = float)
        for pos in list_pos:
            change_pos = pos.copy()
            for i in range(3):
                if self.pos[i] - change_pos[i] > x_cube/2:
                    change_pos[i] += x_cube
                elif self.pos[i] - change_pos[i] < -x_cube/2:
                    change_pos[i] -= x_cube

            delta_r = np.linalg.norm(change_pos - self.pos)
            #print("list_pos = ", list_pos)
            #print("my_pos = ", self.pos)
            #print("delta_r = ", delta_r)
            #if delta_r < 0.5:
            #    print("\n delta_r = ", delta_r, "\n")

            F_mod = 4*epsilon*(-12 * sigma**12 / delta_r**13 + 6 * sigma**6 / delta_r**7)
            F_vec = np.array([0, 0, 0], dtype = float)
            for i in range(3):
                F_vec[i] = -F_mod * (self.pos[i] - change_pos[i]) / delta_r
                F_answ[i] += F_vec[i]

        a = F_answ / mass
        #print('a = ',a)
        #return a
        self.a = a

    def change_pos(self):#, a):
        pos_before = self.pos
        self.pos = 2*self.pos - self.pos_before + self.a * dt**2
        self.pos_before = pos_before
        #pos_new = 2*self.pos - self.pos_before + self.a_new * dt**2
        #self.speed = (self.pos - self.pos_before) / dt
        #self.speed_before = self.speed
        self.speed = (self.pos - self.pos_before) / dt + self.a * dt/2

        #if np.linalg.norm(self.speed) > 10:
        #    print("\n", "speed = ", self.speed, "\n")

        for i in range(3):
            if self.pos[i] > x_cube/2:
                self.pos[i] -= x_cube
                self.pos_before[i] -= x_cube

            if self.pos[i] < -x_cube/2:
                self.pos[i] += x_cube
                self.pos_before[i] += x_cube

    def count_energy(self, list_pos):
        E_kin = mass * (np.linalg.norm(self.speed)**2)/2 #+ np.linalg.norm(self.speed))/2# (self.speed[0]**2 + self.speed[1]**2 + self.speed[2]**2) / 2
        E_pot = 0
        for pos in list_pos:
            change_pos = pos.copy()
            for i in range(3):
                if self.pos[i] - change_pos[i] > x_cube / 2:
                    change_pos[i] += x_cube
                elif self.pos[i] - change_pos[i] < -x_cube / 2:
                    change_pos[i] -= x_cube
            delta_r = np.linalg.norm(change_pos - self.pos)
            E_pot += (4*epsilon * ((sigma / delta_r)**12 - (sigma / delta_r)**6))/2
        return np.array([E_kin, E_pot, E_pot + E_kin], dtype=float) #E_pot + E_kin


class Manager():
    def __init__(self, N):
        self.N = N
        self.t = 0
        self.E = 0

    def random_pos_atom(self):
        V_ball = V_sum_ball/self.N
        list_pos = []
        list_atom = []
        for i in range(self.N):
            main_pos = np.array([0, 0, 0], dtype = float)
            main_speed = np.array([0, 0, 0], dtype = float)
            finish = False
            while not(finish):
                if list_pos == []:
                    finish = True
                for k in range(3):
                    main_pos[k] = random.uniform(-x_cube/2, x_cube/2)
                    main_speed[k] = random.uniform(-speed_max, speed_max)
                for pos in list_pos:
                    pos_c = pos.copy()
                    for i in range(3):
                        if main_pos[i] - pos_c[i] > x_cube / 2:
                            pos_c[i] += x_cube
                        elif main_pos[i] - pos_c[i] < -x_cube / 2:
                            pos_c[i] -= x_cube

                    delta_r = np.linalg.norm(pos_c - main_pos)
                    if delta_r**3 > V_ball:
                        finish = True

            list_pos.append(main_pos)
            list_atom.append(atom(main_pos, main_speed))
        return list_atom

    def main(self, if_E_calculate = False):
        list_atom = self.random_pos_atom()


        if if_E_calculate:
            file_E_kin = open('E_kin.txt', 'w')
            file_E_kin.close()
            file_E_pot = open('E_pot.txt', 'w')
            file_E_pot.close()
            file_E_sum = open('E_kin+E_pot.txt', 'w')
            file_E_sum.close()
            file_V = open('V.txt', 'w')
            file_V.close()
            file_V = open('V.txt', 'a')

        for n in range(N_steps):
            if n % 3 == 0:
                if_E_n = True
            else:
                if_E_n = False

            if n == N_steps - 1:
                print('asdfa')
                xyz = open('mol.xyz', 'w')
                xyz.write(str(N) + '\n')
                xyz.write('Atoms.')
                xyz.write('i = {} \n'.format(n))
                for atom2 in list_atom:
                    line = 'C {} {} {}\n'.format(atom2.pos[0], atom2.pos[1], atom2.pos[2])
                    xyz.write(line)
                xyz.close()

            for i in range(N):
                """
                print("Шаг = ",n)
                print(i, "pos = ", list_atom[i].pos)
                print(i, "speed = ", list_atom[i].speed)
                print(i, "a = ", list_atom[i].a)
                print("\n")"""
                for k in range(3):
                    if abs(list_atom[i].pos[k]) > x_cube/2:
                        print("WARNING WARNING WARNING WARNING WARNING WARNING WARNING WARNING WARNING WARNING WARNING WARNING WARNING WARNING")
                        #print(i, "atom = ", list_atom[i].pos)
                        pass

            #a_list = []
            E = np.array([0, 0, 0], dtype = float)
            list_pos = []
            for atom in list_atom:
                list_pos.append(atom.pos)

            for main_atom in list_atom:
                list_pos_my = list_pos.copy()
                list_pos_my.pop(list_atom.index(main_atom))
                main_atom.count_a(list_pos_my)

                if if_E_calculate and if_E_n:
                    E += main_atom.count_energy(list_pos_my)
            k = 0
            if if_E_calculate and if_E_n:
                file_E_kin = open('E_kin.txt', 'a')
                file_E_pot = open('E_pot.txt', 'a')
                file_E_sum = open('E_kin+E_pot.txt', 'a')
                print(E[0], file=file_E_kin)
                print(E[1], file=file_E_pot)
                print(E[2], file=file_E_sum)

            if n == N_steps - 1 or n == (N_steps*9)//10 or n == (N_steps*8)//10 or n == (N_steps*7)//10:
                print(";", file=file_V)
                for atom in list_atom:
                    speed = np.linalg.norm(atom.speed)
                    print(speed, file=file_V)
                    k += 1

            for atom in list_atom:
                atom.change_pos()


            if n % 10 == 0:
                print(n)



man1 = Manager(N)
man1.main(if_E_calculate = True)

