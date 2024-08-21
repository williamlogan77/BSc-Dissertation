from data_prep import data
import numpy as np
import matplotlib.pyplot as plt

temp, density = "3e8", "1e3"

data_dict = data.get_prepared_data(f"../ppn/nuppn/frames/ppn/i_process_runs/t{temp}_p{density}/x-time.dat")
log_time = np.log10((data_dict["time"][1:] * (24 * 3600 * 365)))
cycle = data_dict["cycle"][1:]
# print(np.log10(data_dict["time"]))
k_B = 1.38064852e-23
T = float(temp)
m_n = 1.674927471e-27

# print(data_dict['']['abundance'][1:])

# fig = plt.figure(figsize=(10, 10))
# neut_density = np.log10(data_dict['NEUT']['abundance']["abundance"][1:] * 6.022e23 * 1e3)
# h_density = np.log10(data_dict['PROT']["abundance"]["abundance"][1:] * 6.022e23 * 1e3)
#
# plt.plot(data_dict["cycle"][1:], neut_density, color="red", label="NEUT")
# # plt.plot(log_time, np.log10(data_dict['C 13']['abundance'][1:]), color="blue", label = "C 13")
# # plt.plot(log_time, np.log10(data_dict['N 13']['abundance'][1:]), color="green", linestyle=":", label ="N 13")
# # plt.plot(log_time, np.log10(data_dict['FE 56']['abundance'][1:]), color="cyan", label="Fe56")
# plt.xticks(np.arange(0, len(data_dict["cycle"]), 20))
# plt.xlabel("Cycle Number")
# plt.ylabel("Log10 Abundance")
# plt.title("Cycle Number vs Neutron Abundance")
# plt.tight_layout()
# plt.legend()
#
# plt.show()

fig, ax = plt.subplots(figsize=(10, 10))
pm154 = data_dict["PM154"]["abundance"][1:]
sm153 = data_dict["SM153"]["abundance"][1:]
sm154 = data_dict["SM154"]["abundance"][1:]
pm153 = data_dict["PM153"]["abundance"][1:]
sm152 = data_dict["SM152"]["abundance"][1:]
sm155 = data_dict["SM155"]["abundance"][1:]

ax.plot(cycle, data_dict["NEUT"]["abundance"]["abundance"][1:], color="red", linestyle=":", label="NEUTRON")
ax.plot(cycle, pm154, color="red", label="PM 154")
ax.plot(cycle, sm153, color="blue", label="SM 153")
ax.plot(cycle, sm154, color="cyan", label="SM 154")
ax.plot(cycle, pm153, color="purple", label="PM 153")
ax.plot(cycle, pm153, color="black", label="PM 153")
ax.plot(cycle, sm152, color="green", label="SM 152")
ax.plot(cycle, sm155, color="orange", label="SM 155")
plt.xticks(np.arange(0, len(data_dict["cycle"]), 50))
ax.set_xlabel("Log10 Time")
ax.set_ylabel("Log10 Abundance")
ax.semilogy()
ax.set_ylim(1e-20)
ax.legend(loc="upper left")
# plt.title("Time(s) vs Neutron Abundance")


ax2 = ax.twinx()

neut_density = (data_dict['NEUT']['abundance']["abundance"][1:] * 6.022e23 * 1e3)

dt = []
for r in range(len(data_dict["time"])):
    if r != 0:
        dt.append(((data_dict["time"][r] - data_dict["time"][r - 1]) * (24 * 3600 * 365)))

neut_exposure = [0]
for q in range(len(data_dict["time"][1:])):
    if q != 0:
        neut_exposure_prime = neut_density[q] * (np.sqrt((2 * k_B * T) / m_n) * 100) * dt[q] * 10e-27
        neut_exposure.append(neut_exposure_prime + neut_exposure[q - 1])


ax2.plot(cycle, neut_exposure, label="Neutron Exposure")
# plt.tight_layout()

ax2.axvline(x=data_dict["cycle"][np.argmax(sm154)], linestyle="--", color="maroon",
            label=f"Sm154 cycle {data_dict['cycle'][np.argmax(sm154)]}")
ax2.axvline(x=data_dict["cycle"][np.argmax(sm155)], linestyle="--", color="magenta",
            label=f"Sm155 cycle {data_dict['cycle'][np.argmax(sm155)]}")

ax2.legend(loc="upper right")
ax2.semilogy()
ax2.set_ylabel("Neutron Exposure $mbarn^-1$")
plt.tight_layout()
plt.show()
