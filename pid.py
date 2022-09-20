#!/usr/bin/env python3
import matplotlib.pyplot as plt
from copy import copy
import time
import math


class PIDController:

    class PIDState:

        def __init__(self, prev_error: float = 0, prev_integral: float = 0, last_updated_s: float = time.time()) -> None:
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

    def __init__(self, kP : float, kI: float, kD: float, target: float) -> None:
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

    def tune(self, kP : float, kI: float, kD: float) -> None:
        self.clear()
        self._kP = kP
        self._kI = kI
        self._kD = kD

    def __str__(self) -> str:
        return self.dump()

    def __repr__(self) -> str:
        return self.__str__()


if __name__ == '__main__':
    targetMax = 3
    # P = .8
    # I = .01
    # D = .05
    P = 1
    I = 0
    D = 0
    p0 = PIDController(P, I, D, targetMax)
    p1 = PIDController(P, I, D, targetMax)
    p2 = PIDController(P, I, D, targetMax)
    input0 = 6
    input1 = 6
    input2 = 6

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

    for i in range(70):
        x.append(i/10)

        target = targetMax*math.pow(2, -i/10)+3

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

        actuals0.append(input0 + output0)
        actuals1.append(input1 + output1)
        actuals2.append(input2 + output2)

        input0 = input0 + output0
        input1 = input1 + output1 + (.1*i**2)/200
        input2 = input2 + output2 - (.1*i**2)/200

    fig, [ax0, ax1, ax2] = plt.subplots(1, 3, constrained_layout=True, figsize=(12, 4))

    ax0.set_title("Ideal Theory")
    ax0.plot(x, inputs0, 'b', label='Inputs', alpha=.5)
    ax0.plot(x, outputs0, 'ro', label='Outputs', alpha=.5)
    ax0.plot(x, actuals0, 'r', label='Actuals')
    ax0.plot(x, targets, 'g--', label='Target')
    ax0.plot(x, [0]*len(x), 'k--')
    ax0.grid()

    ax1.set_title('Correct for "Faster Than Expected"')
    ax1.plot(x, inputs1, 'b', label='Inputs', alpha=.5)
    ax1.plot(x, outputs1, 'ro', label='Outputs', alpha=.5)
    ax1.plot(x, actuals1, 'r', label='Actuals')
    ax1.plot(x, targets, 'g--', label='Target')
    ax1.plot(x, [0]*len(x), 'k--')
    ax1.grid()

    ax2.set_title('Correct for "Slower Than Expected"')
    ax2.plot(x, inputs2, 'b', label='Inputs', alpha=.5)
    ax2.plot(x, outputs2, 'ro', label='Outputs', alpha=.5)
    ax2.plot(x, actuals2, 'r', label='Actuals')
    ax2.plot(x, targets, 'g--', label='Target')
    ax2.plot(x, [0]*len(x), 'k--')
    ax2.grid()

    # plt.legend('Inputs', 'Outputs', 'Actual', 'Target')
    ax2.legend()
    fig.suptitle(f"PID [{P}, {I}, {D}]")
    plt.show()
