import pynq
from pynq import Overlay
from pynq import allocate
import numpy as np
import time


IN = allocate(shape=(2048*2048,), dtype=np.int32)
KER = allocate(shape=(16*16,), dtype=np.int32)
OUT = allocate(shape=(2048*2048,), dtype=np.int32)

IN[:] = range(2048 * 2048)
KER[:] = range(16 * 16)
OUT[:] = [0] * (2048 * 2048)

IN.flush()
KER.flush()

rails = pynq.get_rails()
recorder = pynq.DataRecorder(rails['5V'].power)

overlay = Overlay("/home/xilinx/overlay/design_1.bit")

overlay.convolution_top_0.write(0x10, A.device_address)
overlay.convolution_top_0.write(0x1C, B.device_address)
overlay.convolution_top_0.write(0x28, C.device_address)

exec_time_arr = []
power_arr = []
energy_arr = []
for i in range(10):
	start = time.time()
	overlay.convolution_top_0.write(0x00, 1)

	with recorder.record(0.01):
		while (not overlay.convolution_top_0.read(0x00) & 0x2):
			pass

	end = time.time()

	C.flush()

	exec_time = end - start
	frame = recorder.frame[recorder.frame['Mark'] == 0]
	power = frame["5V_power"].mean()
	energy = power * exec_time

	exec_time_arr.append(exec_time)
	power_arr.append(power)
	energy_arr.append(energy)


print("exec time: " + str(np.avg(exec_time_arr)) +
		"+-" + str(np.std(exec_time_arr)) + " s")
print("power: " + str(np.avg(power_arr)) + "+-" + str(np.std(power_arr)) +
		" (energy: " + str(np.avg(energy_arr)) +
		"+- " + str(np.std(energy_arr)) + ")")

