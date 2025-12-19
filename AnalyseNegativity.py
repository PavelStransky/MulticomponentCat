import numpy as np
import matplotlib.pyplot as plt

PATH = ''

J2 = 4
OMEGA = 1
R = 50

LAMBDA = -np.sqrt(2.0) / 5
DELTA = 0.5

negativity = []
purity = []

times = np.linspace(0, 400, 1000, endpoint=False)

for time in times:
    try:
        data = np.loadtxt(f'{PATH}Wigner_Rabi(2j={J2}, ω={OMEGA}, R={R}, λ={LAMBDA}, δ={DELTA})_{time:.2f}.txt')

    except Exception as e:
        print(f"Error loading {time:.2f}: {e}")
        continue

    y = data[:,2]
    y = np.nan_to_num(y, nan=0.0, posinf=0.0, neginf=0.0)
    w = sum(y)
    y /= w

    negativity.append(0.5 * (sum(abs(y)) - 1))

    print(time)

negativity = np.array(negativity)

plt.figure(figsize=(15, 10))
plt.plot(times, negativity, label='Negativity')
plt.ylim(0, 1)
plt.show()

np.savetxt(f'{PATH}Negativity_Rabi(2j={J2}, ω={OMEGA}, R={R}, λ={LAMBDA}, δ={DELTA}).txt', np.column_stack((times, negativity)))