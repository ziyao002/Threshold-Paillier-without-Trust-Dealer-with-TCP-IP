import socket
from threading import Thread
import multiprocessing
from multiprocessing import Process
import gmpy2
import time
from gmpy2 import mpz
from queue import Queue

# P1:server, P2:server/client, P3:client
server_port1 = 9000
server_port2 = 8000
client_port2 = 8001
client_port31 = 7000
client_port32 = 7001


def rec_msg(tcp_socket, port_num, q12, q13):
    str_list = []
    while True:
        s_data = tcp_socket.recv(1024).decode('gb2312')
        if len(s_data) != 0:
            for item in s_data:
                str_list.append(item)
                if item == 'x':
                    str_data = "".join(str_list[0:-1])
                    # print("rec_data =", int(str_data), "port_num =", port_num)
                    str_list.clear()
                    if port_num == server_port1:
                        q12.put(int(str_data))
                    elif port_num == server_port2:
                        q13.put(int(str_data))
                    else:
                        raise Exception("Port number error!")


def send_msg(tcp_socket, port_num, flag_send_3_to_1, flag_send_3_to_2, data_3_to_1_queue, data_3_to_2_queue):
    while True:
        if port_num == server_port1 and flag_send_3_to_1.value:
            c_data = str(data_3_to_1_queue.get())+'x'
            tcp_socket.send(c_data.encode("gb2312"))
            # print("send_data =", data_3_to_1_queue.value, "port_num =", port_num)
            flag_send_3_to_1.value = 0
        elif port_num == server_port2 and flag_send_3_to_2.value:
            c_data = str(data_3_to_2_queue.get())+'x'
            tcp_socket.send(c_data.encode("gb2312"))
            # print("send_data =", data_3_to_2_queue.value, "port_num =", port_num)
            flag_send_3_to_2.value = 0
        else:
            pass
    # tcp_socket.close()


def client1(flag_send_3_to_1, flag_send_3_to_2, data_3_to_1_queue, data_3_to_2_queue, q31, q32, connect_flag_31):
    host = socket.gethostname()
    client_socket = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
    client_socket.setsockopt(socket.SOL_SOCKET, socket.SO_REUSEADDR, 1)
    client_socket.setsockopt(socket.IPPROTO_TCP, socket.TCP_NODELAY, True)  # disable Nalge
    client_socket.bind((host, client_port31))
    client_socket.connect((host, server_port1))
    connect_flag_31.value = 1
    print("Connect to P1")

    t_send = Thread(target=send_msg, args=(client_socket, server_port1, flag_send_3_to_1, flag_send_3_to_2, data_3_to_1_queue, data_3_to_2_queue))
    t_rec = Thread(target=rec_msg, args=(client_socket, server_port1, q31, q32))

    t_send.start()
    t_rec.start()

    t_send.join()
    t_rec.join()


def client2(flag_send_3_to_1, flag_send_3_to_2, data_3_to_1_queue, data_3_to_2_queue, q31, q32, connect_flag_32):
    host = socket.gethostname()
    client_socket = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
    client_socket.setsockopt(socket.SOL_SOCKET, socket.SO_REUSEADDR, 1)
    client_socket.setsockopt(socket.IPPROTO_TCP, socket.TCP_NODELAY, True)  # disable Nalge
    client_socket.bind((host, client_port32))
    client_socket.connect((host, server_port2))
    # client_socket.connect(('172.21.176.151', server_port2))
    connect_flag_32.value = 1
    print("Connect to P2")

    t_send = Thread(target=send_msg, args=(client_socket, server_port2, flag_send_3_to_1, flag_send_3_to_2, data_3_to_1_queue, data_3_to_2_queue))
    t_rec = Thread(target=rec_msg, args=(client_socket, server_port2, q31, q32))

    t_send.start()
    t_rec.start()

    t_send.join()
    t_rec.join()


class DPaillier:
    def __init__(self, party_index):
        self.KeyLength = mpz(64)
        self.PartyIndex = party_index
        self.PartyNumber = 3
        self.PP = 0
        self.pi = 0
        self.qi = 0
        self.a1 = 0
        self.aa1 = 0
        self.b1 = 0
        self.bb1 = 0
        self.ppi = 0
        self.qqi = 0
        self.c0 = 0
        self.c1 = 0
        self.c2 = 0
        self.cc0 = 0
        self.cc1 = 0
        self.cc2 = 0
        self.Ni = 0
        self.N = 0
        self.Q = 0
        self.gg = 0

        self.generator = mpz(0)
        self.delta = 6
        self.secure_t = 80
        self.secure_l = 80
        self.secure_K = pow(2, 80)
        self.modulus_KN = mpz(0)
        self.modulus_KKN = mpz(0)
        self.beta = mpz(0)
        self.ri = mpz(0)
        self.ri_delta = mpz(0)
        self.thetai = mpz(0)
        self.pq_sum = mpz(0)
        self.phi = mpz(0)
        self.theta = mpz(0)
        self.fi = mpz(0)
        self.r_vk = mpz(0)
        self.vk = mpz(0)
        self.vki = mpz(0)

    def gen_coprime(self, x):
        while True:
            random_state = gmpy2.random_state(int(time.time()*100000))
            coprime = gmpy2.mpz_random(random_state, x)
            if gmpy2.gcd(coprime, x) == 1:
                return coprime

    def pick_pq(self):
        pq = mpz(2)
        random_state = gmpy2.random_state(int(time.time()*100000))
        if self.PartyIndex == 1:
            while gmpy2.f_mod(pq, 4) != 3 or (pq - pow(2, self.KeyLength - 1) <= 0):
                pq = gmpy2.mpz_random(random_state, pow(2, self.KeyLength))
        else:
            while gmpy2.f_mod(pq, 4) != 0 or (pq - pow(2, self.KeyLength - 1) <= 0):
                pq = gmpy2.mpz_random(random_state, pow(2, self.KeyLength))
        return pq

    def pick_pp(self):
        pp = gmpy2.next_prime(pow(gmpy2.mul(self.PartyNumber, gmpy2.mul(3, pow(2, self.KeyLength - 1))), 2))
        return pp

    def coefficient_generation(self):
        random_state = gmpy2.random_state(int(time.time()*100000))
        self.a1 = gmpy2.mpz_random(random_state, self.PP)
        random_state = gmpy2.random_state(int(time.time()*100000))
        self.aa1 = gmpy2.mpz_random(random_state, self.PP)
        random_state = gmpy2.random_state(int(time.time()*100000))
        self.b1 = gmpy2.mpz_random(random_state, self.PP)
        random_state = gmpy2.random_state(int(time.time()*100000))
        self.bb1 = gmpy2.mpz_random(random_state, self.PP)
        random_state = gmpy2.random_state(int(time.time()*100000))
        self.ppi = gmpy2.mpz_random(random_state, self.PP)
        random_state = gmpy2.random_state(int(time.time()*100000))
        self.qqi = gmpy2.mpz_random(random_state, self.PP)
        self.c0 = mpz(0)
        random_state = gmpy2.random_state(int(time.time()*100000))
        self.c1 = gmpy2.mpz_random(random_state, self.PP)
        random_state = gmpy2.random_state(int(time.time()*100000))
        self.c2 = gmpy2.mpz_random(random_state, self.PP)
        random_state = gmpy2.random_state(int(time.time()*100000))
        self.cc0 = gmpy2.mpz_random(random_state, self.PP)
        random_state = gmpy2.random_state(int(time.time()*100000))
        self.cc1 = gmpy2.mpz_random(random_state, self.PP)
        random_state = gmpy2.random_state(int(time.time()*100000))
        self.cc2 = gmpy2.mpz_random(random_state, self.PP)

    def compute_tuple(self):
        self.PP = self.pick_pp()
        # print("pick PP done")
        self.coefficient_generation()
        # print("coefficient generation done")

        pi1 = gmpy2.f_mod((self.pi + self.a1 * 1), self.PP)
        ppi1 = gmpy2.f_mod((self.ppi + self.aa1 * 1), self.PP)
        qi1 = gmpy2.f_mod((self.qi + self.b1 * 1), self.PP)
        qqi1 = gmpy2.f_mod((self.qqi + self.bb1 * 1), self.PP)
        hi1 = gmpy2.f_mod((self.c0 + self.c1 * 1 + self.c2 * 1 * 1), self.PP)
        hhi1 = gmpy2.f_mod((self.cc0 + self.cc1 * 1 + self.cc2 * 1 * 1), self.PP)

        pi2 = gmpy2.f_mod((self.pi + self.a1 * 2), self.PP)
        ppi2 = gmpy2.f_mod((self.ppi + self.aa1 * 2), self.PP)
        qi2 = gmpy2.f_mod((self.qi + self.b1 * 2), self.PP)
        qqi2 = gmpy2.f_mod((self.qqi + self.bb1 * 2), self.PP)
        hi2 = gmpy2.f_mod((self.c0 + self.c1 * 2 + self.c2 * 2 * 2), self.PP)
        hhi2 = gmpy2.f_mod((self.cc0 + self.cc1 * 2 + self.cc2 * 2 * 2), self.PP)

        pi3 = gmpy2.f_mod((self.pi + self.a1 * 3), self.PP)
        ppi3 = gmpy2.f_mod((self.ppi + self.aa1 * 3), self.PP)
        qi3 = gmpy2.f_mod((self.qi + self.b1 * 3), self.PP)
        qqi3 = gmpy2.f_mod((self.qqi + self.bb1 * 3), self.PP)
        hi3 = gmpy2.f_mod((self.c0 + self.c1 * 3 + self.c2 * 3 * 3), self.PP)
        hhi3 = gmpy2.f_mod((self.cc0 + self.cc1 * 3 + self.cc2 * 3 * 3), self.PP)

        return [[pi1, ppi1, qi1, qqi1, hi1, hhi1], [pi2, ppi2, qi2, qqi2, hi2, hhi2], [pi3, ppi3, qi3, qqi3, hi3, hhi3]]

    def send_pq_tuple(self, pq_tuple, send_party_index, flag_send_3_to_1, flag_send_3_to_2, data_3_to_1_queue, data_3_to_2_queue):
        for ctuple in pq_tuple:
            # print("ctuple = ", ctuple, "send_party_index = ", send_party_index)
            self.send_data(ctuple, send_party_index, flag_send_3_to_1, flag_send_3_to_2, data_3_to_1_queue, data_3_to_2_queue)

    def send_data(self, data, party_send_index, flag_send_3_to_1, flag_send_3_to_2, data_3_to_1_queue, data_3_to_2_queue):
        while True:
            if party_send_index == 1 and flag_send_3_to_1.value == 0:
                data_3_to_1_queue.put(data)
                flag_send_3_to_1.value = 1
                break
            elif party_send_index == 2 and flag_send_3_to_2.value == 0:
                data_3_to_2_queue.put(data)
                flag_send_3_to_2.value = 1
                break
            else:
                pass

    def receive_pq_tuple_list(self, self_pq_tuple_list, q31, q32):
        q31_list = []
        q32_list = []
        while True:
            while True:
                while not q31.empty():
                    q31_list.append(mpz(q31.get()))
                    if q31_list:
                        if q31_list[-1] == mpz(11113333):
                            break
                if q31_list:
                    if q31_list[-1] == mpz(11113333):
                        break
            while True:
                while not q32.empty():
                    q32_list.append(mpz(q32.get()))
                    if q32_list:
                        if q32_list[-1] == mpz(22223333):
                            break
                if q32_list:
                    if q32_list[-1] == mpz(22223333):
                        break
            break
        return [q31_list[0:-1], q32_list[0:-1], self_pq_tuple_list[2]]

    def receive_Ni_list(self, self_Ni, q31, q32):
        q31_list = []
        q32_list = []
        while True:
            while True:
                while not q31.empty():
                    q31_list.append(mpz(q31.get()))
                    if q31_list:
                        if q31_list[-1] == 11113333:
                            break
                if q31_list:
                    if q31_list[-1] == 11113333:
                        break
            while True:
                while not q32.empty():
                    q32_list.append(mpz(q32.get()))
                    if q32_list:
                        if q32_list[-1] == 22223333:
                            break
                if q32_list:
                    if q32_list[-1] == 22223333:
                        break
            break
        return [q31_list[0], q32_list[0], self_Ni]

    def send_pq_tuple_list(self, pq_tuple_list, flag_send_3_to_1, flag_send_3_to_2, data_3_to_1_queue, data_3_to_2_queue):
        self.send_pq_tuple(pq_tuple_list[0], 1, flag_send_3_to_1, flag_send_3_to_2, data_3_to_1_queue, data_3_to_2_queue)
        self.send_data(33331111, 1, flag_send_3_to_1, flag_send_3_to_2, data_3_to_1_queue, data_3_to_2_queue)
        self.send_pq_tuple(pq_tuple_list[1], 2, flag_send_3_to_1, flag_send_3_to_2, data_3_to_1_queue, data_3_to_2_queue)
        self.send_data(33332222, 2, flag_send_3_to_1, flag_send_3_to_2, data_3_to_1_queue, data_3_to_2_queue)
        while True:
            if flag_send_3_to_2.value == 0:
                break

    def share_verification(self, received_pq_tuple_list):
        check_flag = 0
        if check_flag == 1:
            raise Exception("Share verification fails!")

    def N_verification(self, Ni_list):
        check_flag = 0
        if check_flag == 1:
            raise Exception("Share verification fails!")

    def compute_Ni(self, received_pq_tuple_list):
        Ni = gmpy2.f_mod(((received_pq_tuple_list[0][0] + received_pq_tuple_list[1][0] + received_pq_tuple_list[2][
            0]) * (received_pq_tuple_list[0][2] + received_pq_tuple_list[1][2] + received_pq_tuple_list[2][2]) + (
                                     received_pq_tuple_list[0][4] + received_pq_tuple_list[1][4] +
                                     received_pq_tuple_list[2][4])), self.PP)
        # print("Ni = ", Ni)
        return Ni

    def compute_N(self, Ni_list):
        L1 = mpz(int((0 - 2) * (0 - 3) / ((1 - 2) * (1 - 3))))
        L2 = mpz(int((0 - 1) * (0 - 3) / ((2 - 1) * (2 - 3))))
        L3 = mpz(int((0 - 1) * (0 - 2) / ((3 - 1) * (3 - 2))))
        self.N = gmpy2.f_mod(gmpy2.mul(Ni_list[0], L1) + gmpy2.mul(Ni_list[1], L2) + gmpy2.mul(Ni_list[2], L3), self.PP)
        # print("Ni_list = ", Ni_list)
        print("Candidate modulus = ", self.N)

    def send_Ni(self, Ni, flag_send_3_to_1, flag_send_3_to_2, data_3_to_1_queue, data_3_to_2_queue):
        self.send_data(Ni, 1, flag_send_3_to_1, flag_send_3_to_2, data_3_to_1_queue, data_3_to_2_queue)
        self.send_data(33331111, 1, flag_send_3_to_1, flag_send_3_to_2, data_3_to_1_queue, data_3_to_2_queue)
        self.send_data(Ni, 2, flag_send_3_to_1, flag_send_3_to_2, data_3_to_1_queue, data_3_to_2_queue)
        self.send_data(33332222, 2, flag_send_3_to_1, flag_send_3_to_2, data_3_to_1_queue, data_3_to_2_queue)
        while True:
            if flag_send_3_to_2.value == 0:
                break

    def receive_gg(self, q31):
        q31_list = []
        while True:
            while not q31.empty():
                q31_list.append(mpz(q31.get()))
                if q31_list:
                    if q31_list[-1] == 11113333:
                        break
            if q31_list:
                if q31_list[-1] == 11113333:
                    break
        ggt = q31_list[0]
        if gmpy2.jacobi(ggt, self.N) == 1:
            self.gg = ggt
        else:
            raise Exception("gg generation Error!")

    def receive_Q_list(self, q31, q32):
        q31_list = []
        q32_list = []
        while True:
            while not q31.empty():
                q31_list.append(mpz(q31.get()))
                if q31_list:
                    if q31_list[-1] == 11113333:
                        break
            if q31_list:
                if q31_list[-1] == 11113333:
                    break
        while True:
            while not q32.empty():
                q32_list.append(mpz(q32.get()))
                if q32_list:
                    if q32_list[-1] == 22223333:
                        break
            if q32_list:
                if q32_list[-1] == 22223333:
                    break
        return [q31_list[0], q32_list[0], self.Q]

    def biprimality_check(self, flag_send_3_to_1, flag_send_3_to_2, data_3_to_1_queue, data_3_to_2_queue, q31, q32):
        self.receive_gg(q31)
        self.Q = gmpy2.powmod(self.gg, gmpy2.f_div(self.pi + self.qi, 4), self.N)

        self.send_data(self.Q, 1, flag_send_3_to_1, flag_send_3_to_2, data_3_to_1_queue, data_3_to_2_queue)
        self.send_data(33331111, 1, flag_send_3_to_1, flag_send_3_to_2, data_3_to_1_queue, data_3_to_2_queue)
        self.send_data(self.Q, 2, flag_send_3_to_1, flag_send_3_to_2, data_3_to_1_queue, data_3_to_2_queue)
        self.send_data(33332222, 2, flag_send_3_to_1, flag_send_3_to_2, data_3_to_1_queue, data_3_to_2_queue)
        while True:
            if flag_send_3_to_2.value == 0:
                break
        Q_list = self.receive_Q_list(q31, q32)
        # print("Q_list = ", Q_list)
        # print("Q_list = ", Q_list)

        Q1 = Q_list[0]
        Q2 = Q_list[1]
        Q3 = Q_list[2]
        Q2_inv = gmpy2.invert(Q2, self.N)
        Q3_inv = gmpy2.invert(Q3, self.N)

        biprimality_check = gmpy2.f_mod((Q1 * Q2_inv * Q3_inv), self.N) == gmpy2.f_mod(mpz(1), self.N) or gmpy2.f_mod(
            (Q1 * Q2_inv * Q3_inv), self.N) == gmpy2.f_mod(mpz(-1), self.N)
        return biprimality_check

    def send_pq_sum(self, pq_sum, flag_send_3_to_1, flag_send_3_to_2, data_3_to_1_queue, data_3_to_2_queue):
        self.send_data(pq_sum, 1, flag_send_3_to_1, flag_send_3_to_2, data_3_to_1_queue, data_3_to_2_queue)
        self.send_data(33331111, 1, flag_send_3_to_1, flag_send_3_to_2, data_3_to_1_queue, data_3_to_2_queue)
        self.send_data(pq_sum, 2, flag_send_3_to_1, flag_send_3_to_2, data_3_to_1_queue, data_3_to_2_queue)
        self.send_data(33332222, 2, flag_send_3_to_1, flag_send_3_to_2, data_3_to_1_queue, data_3_to_2_queue)
        while True:
            if flag_send_3_to_2.value == 0:
                break

    def receive_pq_sum_list(self, q31, q32):
        q31_list = []
        q32_list = []
        while True:
            while not q31.empty():
                q31_list.append(mpz(q31.get()))
                if q31_list:
                    if q31_list[-1] == 11113333:
                        break
            if q31_list:
                if q31_list[-1] == 11113333:
                    break
        while True:
            while not q32.empty():
                q32_list.append(mpz(q32.get()))
                if q32_list:
                    if q32_list[-1] == 22223333:
                        break
            if q32_list:
                if q32_list[-1] == 22223333:
                    break
        return [q31_list[0], q32_list[0], self.pq_sum]

    def send_thetai(self, thetai, flag_send_3_to_1, flag_send_3_to_2, data_3_to_1_queue, data_3_to_2_queue):
        self.send_data(thetai, 1, flag_send_3_to_1, flag_send_3_to_2, data_3_to_1_queue, data_3_to_2_queue)
        self.send_data(33331111, 1, flag_send_3_to_1, flag_send_3_to_2, data_3_to_1_queue, data_3_to_2_queue)
        self.send_data(thetai, 2, flag_send_3_to_1, flag_send_3_to_2, data_3_to_1_queue, data_3_to_2_queue)
        self.send_data(33332222, 2, flag_send_3_to_1, flag_send_3_to_2, data_3_to_1_queue, data_3_to_2_queue)
        while True:
            if flag_send_3_to_2.value == 0:
                break

    def receive_thetai_list(self, q31, q32):
        q31_list = []
        q32_list = []
        while True:
            while not q31.empty():
                q31_list.append(mpz(q31.get()))
                if q31_list:
                    if q31_list[-1] == 11113333:
                        break
            if q31_list:
                if q31_list[-1] == 11113333:
                    break
        while True:
            while not q32.empty():
                q32_list.append(mpz(q32.get()))
                if q32_list:
                    if q32_list[-1] == 22223333:
                        break
            if q32_list:
                if q32_list[-1] == 22223333:
                    break
        return [q31_list[0], q32_list[0], self.thetai]

    def send_rh(self, r31, r32, h31, h32, flag_send_3_to_1, flag_send_3_to_2, data_3_to_1_queue, data_3_to_2_queue):
        self.send_data(r31, 1, flag_send_3_to_1, flag_send_3_to_2, data_3_to_1_queue, data_3_to_2_queue)
        self.send_data(h31, 1, flag_send_3_to_1, flag_send_3_to_2, data_3_to_1_queue, data_3_to_2_queue)
        self.send_data(33331111, 1, flag_send_3_to_1, flag_send_3_to_2, data_3_to_1_queue, data_3_to_2_queue)
        self.send_data(r32, 2, flag_send_3_to_1, flag_send_3_to_2, data_3_to_1_queue, data_3_to_2_queue)
        self.send_data(h32, 2, flag_send_3_to_1, flag_send_3_to_2, data_3_to_1_queue, data_3_to_2_queue)
        self.send_data(33332222, 2, flag_send_3_to_1, flag_send_3_to_2, data_3_to_1_queue, data_3_to_2_queue)
        while True:
            if flag_send_3_to_2.value == 0:
                break

    def receive_rh_list(self, q31, q32):
        q31_list = []
        q32_list = []
        while True:
            while not q31.empty():
                q31_list.append(mpz(q31.get()))
                if q31_list:
                    if q31_list[-1] == 11113333:
                        break
            if q31_list:
                if q31_list[-1] == 11113333:
                    break
        while True:
            while not q32.empty():
                q32_list.append(mpz(q32.get()))
                if q32_list:
                    if q32_list[-1] == 22223333:
                        break
            if q32_list:
                if q32_list[-1] == 22223333:
                    break
        return [q31_list[0:-1], q32_list[0:-1]]

    def Ri_sharing(self, r, modulus, flag_send_3_to_1, flag_send_3_to_2, data_3_to_1_queue, data_3_to_2_queue, q31, q32):
        random_state = gmpy2.random_state(int(time.time() * 100000))
        a1 = gmpy2.mpz_random(random_state, modulus)
        random_state = gmpy2.random_state(int(time.time() * 100000))
        c1 = gmpy2.mpz_random(random_state, modulus)
        random_state = gmpy2.random_state(int(time.time() * 100000))
        c2 = gmpy2.mpz_random(random_state, modulus)

        r31 = r + a1 * 1
        r32 = r + a1 * 2
        r33 = r + a1 * 3
        h31 = c1 * 1 + c2 * 1 * 1
        h32 = c1 * 2 + c2 * 2 * 2
        h33 = c1 * 3 + c2 * 3 * 3

        self.send_rh(r31, r32, h31, h32, flag_send_3_to_1, flag_send_3_to_2, data_3_to_1_queue, data_3_to_2_queue)
        rh_list = self.receive_rh_list(q31, q32)
        Ri = rh_list[0][0] + rh_list[0][1] +r33 + rh_list[1][0] + rh_list[1][1] + h33

        return Ri

    def distributed_RSA_modulus_generation(self, flag_send_3_to_1, flag_send_3_to_2, data_3_to_1_queue, data_3_to_2_queue, q31, q32):
        while True:
            # print("Distributed RSA modulus generation start")
            self.pi = self.pick_pq()
            self.qi = self.pick_pq()
            # print("pick pq done")
            pq_tuple_list = self.compute_tuple()
            # print("send_pq_tuple_list = ", pq_tuple_list)
            # print("compute pq tuple done")
            self.send_pq_tuple_list(pq_tuple_list, flag_send_3_to_1, flag_send_3_to_2, data_3_to_1_queue, data_3_to_2_queue)
            # print("send pq tuple list done")
            received_pq_tuple_list = self.receive_pq_tuple_list(pq_tuple_list, q31, q32)
            # print("received_pq_tuple_list", received_pq_tuple_list)
            # print("receive pq tuple list done")
            self.share_verification(received_pq_tuple_list)
            Ni = self.compute_Ni(received_pq_tuple_list)
            # print("compute Ni done")
            self.send_Ni(Ni, flag_send_3_to_1, flag_send_3_to_2, data_3_to_1_queue, data_3_to_2_queue)
            # print("send Ni done")
            Ni_list = self.receive_Ni_list(Ni, q31, q32)
            # print("receive Ni list done")
            # print("Ni_list = ", Ni_list)
            self.N_verification(Ni_list)
            self.compute_N(Ni_list)
            # print("compute N done")
            if self.biprimality_check(flag_send_3_to_1, flag_send_3_to_2, data_3_to_1_queue, data_3_to_2_queue, q31, q32):
                break
            # print("biprimality check done")
            # self.start_sync(q31, q32)

    def distributed_Paillier_key_generation(self, flag_send_3_to_1, flag_send_3_to_2, data_3_to_1_queue, data_3_to_2_queue, q31, q32):
        self.generator = self.N + 1
        self.modulus_KN = self.secure_K * self.N
        self.modulus_KKN = self.secure_K * self.secure_K * self.N

        random_state = gmpy2.random_state(int(time.time() * 100000))
        self.beta = gmpy2.mpz_random(random_state, self.modulus_KN)
        random_state = gmpy2.random_state(int(time.time() * 100000))
        self.ri = gmpy2.mpz_random(random_state, self.modulus_KKN)
        self.ri_delta = self.delta * self.ri

        self.pq_sum = self.pi + self.qi
        self.send_pq_sum(self.pq_sum, flag_send_3_to_1, flag_send_3_to_2, data_3_to_1_queue, data_3_to_2_queue)
        pq_sum_list = self.receive_pq_sum_list(q31, q32)
        self.phi = self.N + 1 - pq_sum_list[0] - pq_sum_list[1] - pq_sum_list[2]

        self.thetai = self.delta * self.phi * self.beta + self.N * self.delta * self.ri
        self.send_thetai(self.thetai, flag_send_3_to_1, flag_send_3_to_2, data_3_to_1_queue, data_3_to_2_queue)
        thetai_list = self.receive_thetai_list(q31, q32)
        self.theta = thetai_list[0] + thetai_list[1] + thetai_list[2]

        self.fi = self.N * self.Ri_sharing(self.ri_delta, self.modulus_KKN, flag_send_3_to_1, flag_send_3_to_2, data_3_to_1_queue, data_3_to_2_queue, q31, q32) - self.theta

        # verification key
        self.r_vk = self.gen_coprime(self.N * self.N)
        self.vk = gmpy2.powmod(self.r_vk, 2, self.N * self.N)
        self.vki = gmpy2.powmod(self.vk, self.delta * self.fi, self.N * self.N)

    def encrypt(self, massage):
        random = self.gen_coprime(self.N)
        ciphertext = gmpy2.f_mod(gmpy2.mul(gmpy2.powmod(self.generator, massage, self.N * self.N), gmpy2.powmod(random, self.N, self.N * self.N)), self.N * self.N)
        return ciphertext

    def partial_decrypt(self, ciphertext):
        partial_decryption = gmpy2.powmod(ciphertext, 2 * self.delta * self.fi, self.N * self.N)
        return partial_decryption

    def L_function(self, u):
        return (u-1)//self.N

    def combine_partial_decrypt(self, partial_decryption_list):
        lamda1 = mpz(3)
        lamda2 = mpz(-3)
        lamda3 = mpz(1)

        product1 = gmpy2.powmod(partial_decryption_list[0], 2 * self.delta * lamda1, self.N * self.N)
        product2 = gmpy2.powmod(partial_decryption_list[1], 2 * self.delta * lamda2, self.N * self.N)
        product3 = gmpy2.powmod(partial_decryption_list[2], 2 * self.delta * lamda3, self.N * self.N)
        product = gmpy2.f_mod(gmpy2.mul(gmpy2.mul(product1, product2), product3), self.N * self.N)

        Inv_temp = gmpy2.invert(gmpy2.f_mod(-4 * self.delta * self.delta * self.theta, self.N), self.N)
        M = gmpy2.f_mod(self.L_function(product) * Inv_temp, self.N)

        return M

    def send_ciphertext(self, ciphertext, flag_send_3_to_1, flag_send_3_to_2, data_3_to_1_queue, data_3_to_2_queue):
        self.send_data(ciphertext, 1, flag_send_3_to_1, flag_send_3_to_2, data_3_to_1_queue, data_3_to_2_queue)
        self.send_data(22221111, 1, flag_send_3_to_1, flag_send_3_to_2, data_3_to_1_queue, data_3_to_2_queue)
        self.send_data(ciphertext, 3, flag_send_3_to_1, flag_send_3_to_2, data_3_to_1_queue, data_3_to_2_queue)
        self.send_data(22223333, 3, flag_send_3_to_1, flag_send_3_to_2, data_3_to_1_queue, data_3_to_2_queue)
        while True:
            if flag_send_3_to_2.value == 0:
                break

    def send_partial_decryption(self, partial_decryption, flag_send_3_to_1, flag_send_3_to_2, data_3_to_1_queue, data_3_to_2_queue):
        self.send_data(partial_decryption, 1, flag_send_3_to_1, flag_send_3_to_2, data_3_to_1_queue, data_3_to_2_queue)
        self.send_data(33331111, 1, flag_send_3_to_1, flag_send_3_to_2, data_3_to_1_queue, data_3_to_2_queue)
        self.send_data(partial_decryption, 2, flag_send_3_to_1, flag_send_3_to_2, data_3_to_1_queue, data_3_to_2_queue)
        self.send_data(33332222, 2, flag_send_3_to_1, flag_send_3_to_2, data_3_to_1_queue, data_3_to_2_queue)
        while True:
            if flag_send_3_to_2.value == 0:
                break

    def rec_partial_decryption_list(self, partial_decryption, q31, q32):
        q31_list = []
        q32_list = []
        while True:
            while not q31.empty():
                q31_list.append(mpz(q31.get()))
                if q31_list:
                    if q31_list[-1] == 11113333:
                        break
            if q31_list:
                if q31_list[-1] == 11113333:
                    break
        while True:
            while not q32.empty():
                q32_list.append(mpz(q32.get()))
                if q32_list:
                    if q32_list[-1] == 22223333:
                        break
            if q32_list:
                if q32_list[-1] == 22223333:
                    break
        return [q31_list[0], q32_list[0], partial_decryption]

    def rec_ciphertext(self, q31):
        q31_list = []
        while True:
            while not q31.empty():
                q31_list.append(mpz(q31.get()))
                if q31_list:
                    if q31_list[-1] == 11113333:
                        break
            if q31_list:
                if q31_list[-1] == 11113333:
                    break
        ciphertext = q31_list[0]
        return ciphertext


if __name__ == "__main__":
    flag_send_3_to_1 = multiprocessing.Value('l', 0)
    flag_send_3_to_2 = multiprocessing.Value('l', 0)
    data_3_to_1_queue = multiprocessing.Queue()
    data_3_to_2_queue = multiprocessing.Queue()
    q31 = multiprocessing.Queue()
    q32 = multiprocessing.Queue()
    connect_flag_31 = multiprocessing.Value('h', 0)
    connect_flag_32 = multiprocessing.Value('h', 0)

    client1_process = Process(target=client1, args=(flag_send_3_to_1, flag_send_3_to_2, data_3_to_1_queue, data_3_to_2_queue, q31, q32, connect_flag_31))
    client2_process = Process(target=client2, args=(flag_send_3_to_1, flag_send_3_to_2, data_3_to_1_queue, data_3_to_2_queue, q31, q32, connect_flag_32))

    client1_process.start()
    client2_process.start()

    # wait for connection
    while True:
        if connect_flag_31.value and connect_flag_32.value:
            break

    # distributed Paillier key generation
    start = time.time()
    d_paillier = DPaillier(3)
    d_paillier.distributed_RSA_modulus_generation(flag_send_3_to_1, flag_send_3_to_2, data_3_to_1_queue, data_3_to_2_queue, q31, q32)
    stop = time.time()
    print("RSA modulus generation success")
    print("modulus = ", d_paillier.N)
    print("duration = ", stop - start, "seconds")
    d_paillier.distributed_Paillier_key_generation(flag_send_3_to_1, flag_send_3_to_2, data_3_to_1_queue, data_3_to_2_queue, q31, q32)

    # decryption and encryption
    ciphertext = d_paillier.rec_ciphertext(q31)
    print("ciphertext = ", ciphertext)
    partial_decryption = d_paillier.partial_decrypt(ciphertext)
    print("partial decryption = ", partial_decryption)
    d_paillier.send_partial_decryption(partial_decryption, flag_send_3_to_1, flag_send_3_to_2, data_3_to_1_queue, data_3_to_2_queue)
    partial_decryption_list = d_paillier.rec_partial_decryption_list(partial_decryption, q31, q32)
    print("partial_decryption_list = ", partial_decryption_list)
    decrypted_massage = d_paillier.combine_partial_decrypt(partial_decryption_list)
    print("decrypted massage = ", decrypted_massage)


