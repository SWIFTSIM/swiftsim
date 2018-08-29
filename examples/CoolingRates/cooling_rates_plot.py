import matplotlib.pyplot as plt

elements = 11
u = []
cooling_rate = [[] for i in range(elements+1)]
file_in = open('cooling_output.dat', 'r')
for line in file_in:
	data = line.split()
	u.append(float(data[0]))
	cooling_rate[0].append(-float(data[1]))

file_in.close()

for elem in range(elements):
	file_in = open('cooling_element_'+str(elem)+'.dat','r')
	for line in file_in:
        	data = line.split()
        	cooling_rate[elem+1].append(-float(data[0]))
	file_in.close()

ax = plt.subplot(111)
p0, = plt.loglog(u, cooling_rate[0], linewidth = 0.5, color = 'k', label = 'Total')
p1, = plt.loglog(u, cooling_rate[1], linewidth = 0.5, color = 'k', linestyle = '--', label = 'H + He')
p2, = plt.loglog(u, cooling_rate[3], linewidth = 0.5, color = 'b', label = 'Carbon')
p3, = plt.loglog(u, cooling_rate[4], linewidth = 0.5, color = 'g', label = 'Nitrogen')
p4, = plt.loglog(u, cooling_rate[5], linewidth = 0.5, color = 'r', label = 'Oxygen')
p5, = plt.loglog(u, cooling_rate[6], linewidth = 0.5, color = 'c', label = 'Neon')
p6, = plt.loglog(u, cooling_rate[7], linewidth = 0.5, color = 'm', label = 'Magnesium')
p7, = plt.loglog(u, cooling_rate[8], linewidth = 0.5, color = 'y', label = 'Silicon')
p8, = plt.loglog(u, cooling_rate[9], linewidth = 0.5, color = 'lightgray', label = 'Sulphur')
p9, = plt.loglog(u, cooling_rate[10], linewidth = 0.5, color = 'olive', label = 'Calcium')
p10, = plt.loglog(u, cooling_rate[11], linewidth = 0.5, color = 'saddlebrown', label = 'Iron')
ax.set_position([0.15,0.15,0.75,0.75])
plt.xlim([1e12,1e17])
plt.ylim([1e-24,1e-21])
plt.xlabel("Internal energy ${\\rm{[erg \cdot g^{-1}]}}$", fontsize = 14)
plt.ylabel("${\Lambda/n_h^2 }$ ${\\rm{[erg \cdot cm^3 \cdot s^{-1}]}}$", fontsize = 14)
plt.legend(handles = [p0,p1,p2,p3,p4,p5,p6,p7,p8,p9,p10])
plt.show()
