#!/usr/bin/env python3
from random import random

import matplotlib.pyplot as plt
from copy import copy
import time
import math

import numpy as np


class PIDController:
    class PIDState:

        def __init__(self, prev_error: float = 0, prev_integral: float = 0,
                     last_updated_s: float = time.time()) -> None:
            self._prev_error = prev_error
            self._prev_integral = prev_integral
            self._last_updated_s = last_updated_s
            self._output = 0

        def dump(self) -> str:
            return f"output = {self._output}, prev_error = {self._prev_error}, prev_integral = {self._prev_integral}, last_updated_s = {self._last_updated_s}"

        def clear(self) -> None:
            self._output = 0
            self._prev_error = 0
            self._prev_integral = 0
            self._last_updated_s = time.time()

        def __str__(self) -> str:
            return self.dump()

        def __repr__(self) -> str:
            return self.__str__()

    def __init__(self, kP: float, kI: float, kD: float, target: float) -> None:
        self._kP = kP
        self._kI = kI
        self._kD = kD
        self._target = target
        self._state = PIDController.PIDState()

    def dump(self) -> str:
        return f"kP = {self._kP}, kI = {self._kI}, kD = {self._kD}, target = {self._target}, {self._state}"

    def update(self, current: float, dt: float = 0) -> float:
        state = self._state

        now_ns = time.time()
        # not used, use the value passed from client.
        _dt = now_ns - state._last_updated_s

        # if not math.isclose(self._target, current, abs_tol=0.5):
        error = self._target - current
        # else:
        #     error = 0

        porportional = self._kP * error

        # if error == 0:
        #     integral = 0
        # else:
        integral = state._prev_integral + self._kI * (error * (1 if dt == 0 else dt))

        derivative = self._kD * ((error - state._prev_error) / (1 if dt == 0 else dt))

        output = porportional + integral + derivative

        state._prev_error = error
        state._prev_integral = integral
        state._last_updated_s = time.time()
        state._output = output

        return output

    def get(self) -> PIDState:
        return copy(self._state)

    def clear(self) -> None:
        self._state.clear()

    def tune(self, kP: float, kI: float, kD: float) -> None:
        self.clear()
        self._kP = kP
        self._kI = kI
        self._kD = kD

    def __str__(self) -> str:
        return self.dump()

    def __repr__(self) -> str:
        return self.__str__()


def update_max_error():
    global max_err_0, max_err_0x, max_err_0y, max_err_1, max_err_1x, max_err_1y, max_err_2, max_err_2x, max_err_2y
    if err_0 > max_err_0:
        max_err_0 = err_0
        max_err_0x = i
        max_err_0y = actuals0[-1]
    if err_1 > max_err_1:
        max_err_1 = err_1
        max_err_1x = i
        max_err_1y = actuals1[-1]
    if err_2 > max_err_2:
        max_err_2 = err_2
        max_err_2x = i
        max_err_2y = actuals2[-1]


def init():
    global p0, p1, p2, input0, input1, input2, output0, output1, output2, inputs0, outputs0, targets, actuals0, x, inputs1, outputs1, actuals1, inputs2, outputs2, actuals2, max_err_0, max_err_0x, max_err_0y, max_err_1, max_err_1x, max_err_1y, max_err_2, max_err_2x, max_err_2y
    p0 = PIDController(P, I, D, 0)
    p1 = PIDController(P, I, D, 0)
    p2 = PIDController(P, I, D, 0)
    input0 = 0
    input1 = 0
    input2 = 0
    output0 = p0.update(input0)
    output1 = p1.update(input1)
    output2 = p2.update(input2)
    inputs0 = []
    outputs0 = []
    targets = []
    actuals0 = []
    x = []
    inputs1 = []
    outputs1 = []
    actuals1 = []
    inputs2 = []
    outputs2 = []
    actuals2 = []
    max_err_0 = 0
    max_err_0x = 0
    max_err_0y = 0
    max_err_1 = 0
    max_err_1x = 0
    max_err_1y = 0
    max_err_2 = 0
    max_err_2x = 0
    max_err_2y = 0


def plotter():
    max_abs_error = max([max_err_0, max_err_1, max_err_2])*1.1

    fig, axs = plt.subplots(2, 3, constrained_layout=True, figsize=(18, 10))
    axs[0, 0].set_title("Ideal With Noise")
    # axs[0, 0].plot(x, inputs0, 'b', label='Inputs', alpha=.5)
    # ax0.plot(x, outputs0, 'ro', label='Outputs', alpha=.5)
    axs[0, 0].plot(x, actuals0, 'r', label='Actuals')
    axs[0, 0].plot(x, targets, 'g--', label='Target')
    axs[0, 0].plot(x, [0] * len(x), 'k--')
    axs[0, 0].plot(max_err_0x, max_err_0y, 'kx')
    axs[0, 0].grid()

    axs[1, 0].set_title(f"Error (max {round(max_err_0)})")
    axs[1, 0].plot(x, [a - t for a, t in zip(actuals0, targets)], 'b', label='Error', alpha=.5)
    # ax0.plot(x, outputs0, 'ro', label='Outputs', alpha=.5)
    axs[1, 0].plot(x, [0] * len(x), 'k--')
    axs[1, 0].plot([max_err_0x]*10, np.linspace(-max_err_0, max_err_0, 10), 'b--')
    axs[1, 0].set_ylim([-max_abs_error, max_abs_error])
    axs[1, 0].grid()

    axs[0, 1].set_title('"Faster Than Expected" with Noise')
    # axs[0,1].plot(x, inputs1, 'b', label='Inputs', alpha=.5)
    # ax1.plot(x, outputs1, 'ro', label='Outputs', alpha=.5)
    axs[0, 1].plot(x, actuals1, 'r', label='Actuals')
    axs[0, 1].plot(x, targets, 'g--', label='Target')
    axs[0, 1].plot(x, [0] * len(x), 'k--')
    axs[0, 1].plot(max_err_1x, max_err_1y, 'kx')
    axs[0, 1].grid()

    axs[1, 1].set_title(f"Error (max {round(max_err_1)})")
    axs[1, 1].plot(x, [a - t for a, t in zip(actuals1, targets)], 'b', label='Error', alpha=.5)
    # ax0.plot(x, outputs0, 'ro', label='Outputs', alpha=.5)
    axs[1, 1].plot(x, [0] * len(x), 'k--')
    axs[1, 1].plot([max_err_1x]*10, np.linspace(-max_err_1, max_err_1, 10), 'b--')
    axs[1, 1].set_ylim([-max_abs_error, max_abs_error])
    axs[1, 1].grid()

    axs[0, 2].set_title('"Slower Than Expected" with Noise')
    # axs[0,2].plot(x, inputs2, 'b', label='Inputs', alpha=.5)
    # ax2.plot(x, outputs2, 'ro', label='Outputs', alpha=.5)
    axs[0, 2].plot(x, actuals2, 'r', label='Actuals')
    axs[0, 2].plot(x, targets, 'g--', label='Target')
    axs[0, 2].plot(x, [0] * len(x), 'k--')
    axs[0, 2].plot(max_err_2x, max_err_2y, 'kx')
    axs[0, 2].grid()

    axs[1, 2].set_title(f"Error (max {round(max_err_2)})")
    axs[1, 2].plot(x, [a - t for a, t in zip(actuals2, targets)], 'b', label='Error', alpha=.5)
    # ax0.plot(x, outputs0, 'ro', label='Outputs', alpha=.5)
    axs[1, 2].plot(x, [0] * len(x), 'k--')
    axs[1, 2].plot([max_err_2x]*10, np.linspace(-max_err_2, max_err_2, 10), 'b--')
    axs[1, 2].set_ylim([-max_abs_error, max_abs_error])
    axs[1, 2].grid()
    # plt.legend('Inputs', 'Outputs', 'Actual', 'Target')
    axs[0, 2].legend()
    axs[1, 2].legend()
    fig.suptitle(f"PID [{P}, {I}, {D}]")
    plt.show()


def insert_noise():
    global input0, input1, input2
    input0 = actuals0[-1] + 2 * f * (random() - .5)
    # input1 = input1 + output1 + (.01*i**2)/20000
    input1 = actuals1[-1] + f * random()
    # input2 = input2 + output2 - (.01*i**2)/20000
    input2 = actuals2[-1] - f * random()


def update_target_and_reaction():
    global output0, output1, output2, err_0, err_1, err_2
    x.append(i)
    target = max(0, -2066.8 + 1.9643 * i \
                 - 3.9390E-06 * math.pow(i, 2) \
                 + 1.0293E-11 * math.pow(i, 3) \
                 - 2.1983E-17 * math.pow(i, 4) \
                 + 2.8266E-23 * math.pow(i, 5) \
                 - 1.5692E-29 * math.pow(i, 6))
    p0._target = target
    p1._target = target
    p2._target = target
    targets.append(target)
    inputs0.append(input0)
    inputs1.append(input1)
    inputs2.append(input2)
    output0 = p0.update(input0)
    output1 = p1.update(input1)
    output2 = p1.update(input2)
    outputs0.append(output0)
    outputs1.append(output1)
    outputs2.append(output2)
    actuals0.append(max(0, input0 + output0))
    actuals1.append(max(0, input1 + output1))
    actuals2.append(max(0, input2 + output2))
    err_0 = abs(actuals0[-1] - targets[-1])
    err_1 = abs(actuals1[-1] - targets[-1])
    err_2 = abs(actuals2[-1] - targets[-1])


if __name__ == '__main__':
    # P = .8
    # I = .01
    # D = .05
    P = 0.025
    I = 0.001
    D = 0.0001

    init()

    max_clicks = int(4.6 * math.pow(10, 5))
    for i in range(max_clicks):
        update_target_and_reaction()

        update_max_error()

        f = 40
        insert_noise()

    print(f"MAX_ERR_0 = {max_err_0}")
    print(f"MAX_ERR_1 = {max_err_1}")
    print(f"MAX_ERR_2 = {max_err_2}")

    plotter()

