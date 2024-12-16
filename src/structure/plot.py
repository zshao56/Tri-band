# plot the txt file
import matplotlib.pyplot as plt
    
with open('results_N2.txt', 'r') as f:
    lines = f.readlines()[1:]  # Skip the first line
        
    wvl = []
    y = []
    for line in lines:
        data = line.split()
        wvl.append(float(data[0]))
        .append(float(data[1]))
            
    plt.plot(x, y)
    plt.show()