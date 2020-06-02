import matplotlib.pyplot as plt
import numpy as np

_direction = [[0, 1], [1, 0], [0, 1]]
_eps = 0.01


def foo(x, xd, y, yd, nd, lam=0, dx=0, operation=1):
    n = 8
    px = x + (lam - operation * dx) * xd / nd
    py = y + (lam - operation * dx) * yd / nd
    return 5 * (px - n) ** 2 + px * py + 3 * py ** 2


def setup(start, direction, sb=True):
    x, y = start[0], start[1]
    xd, yd = direction[0], direction[1]
    nd = (xd ** 2 + yd ** 2) ** 0.5
    ns = (x ** 2 + y ** 2) ** 0.5
    dx = 0.1 * (ns / nd)
    if sb:
        return x, y, xd, yd, nd, ns, dx
    else:
        return x, y, xd, yd, nd


def sven(start, direction, lam=0, dsk_=False):
    x0, y0, xd, yd, normalized_direction, normalized_start, dx = setup(start, direction)

    values = [foo(x0, xd, y0, yd, normalized_direction, lam)]
    lambdas = [lam]
    if foo(x0, xd, y0, yd, normalized_direction, lam, dx) < \
            foo(x0, xd, y0, yd, normalized_direction, lam) < \
            foo(x0, xd, y0, yd, normalized_direction, lam, dx, -1):
        values.append(foo(x0, xd, y0, yd, normalized_direction, lam, dx))
        lambdas.append(lam - dx)
    elif foo(x0, xd, y0, yd, normalized_direction, lam, dx) > \
            foo(x0, xd, y0, yd, normalized_direction, lam) > \
            foo(x0, xd, y0, yd, normalized_direction, lam, dx, -1):
        values.append(foo(x0, xd, y0, yd, normalized_direction, lam, dx, -1))
        lambdas.append(lam + dx)
    elif foo(x0, xd, y0, yd, normalized_direction, lam, dx) > \
            foo(x0, xd, y0, yd, normalized_direction, lam, ) < \
            foo(x0, xd, y0, yd, normalized_direction, lam, dx, -1):
        if not dsk_:
            return lam - dx, lam + dx
        else:
            return lam-dx, lam, lam+dx
    operation = 1 if lambdas[1] > lambdas[0] else -1
    i = 1
    while values[i] < values[i - 1]:
        step = lambdas[i] + operation * dx * (2 ** i)
        lambdas.append(step)

        values.append(foo(x0, xd, y0, yd, normalized_direction, step, -operation))
        i += 1

    last = (lambdas[i], (lambdas[i] + lambdas[i - 1])/2, lambdas[i - 1], lambdas[i - 2])
    eval_lambdas = []
    for i in last:
        eval_lambdas.append(foo(x0, xd, y0, yd, normalized_direction, i))

    minimal = eval_lambdas.index(min(eval_lambdas))

    if not dsk_:
        return sorted([last[minimal - 1], last[minimal + 1]])
    else:
        return sorted((last[minimal - 1], last[minimal], last[minimal + 1]))


def dichotomy(sven_, start, direction):
    x0, y0, xd, yd, nd = setup(start, direction, sb=False)
    long = sven_[1] - sven_[0]
    middle = (sven_[1] + sven_[0]) / 2
    left = (sven_[0] + middle) / 2
    right = (middle + sven_[1]) / 2
    while long > _eps:
        if foo(x0, xd, y0, yd, nd, left) < \
                foo(x0, xd, y0, yd, nd, middle) < \
                foo(x0, xd, y0, yd, nd, right):
            right = (left + middle) / 2
            middle = left
            left = (sven_[0] + middle) / 2
        elif foo(x0, xd, y0, yd, nd, left) > \
                foo(x0, xd, y0, yd, nd, middle) > \
                foo(x0, xd, y0, yd, nd, right):
            left = (middle + right) / 2
            middle = right
            right = (middle + sven_[1]) / 2
        elif foo(x0, xd, y0, yd, nd, left) > \
                foo(x0, xd, y0, yd, nd, middle) < \
                foo(x0, xd, y0, yd, nd, right):
            left = (left + middle) / 2
            middle = middle
            right = (right + middle) / 2
        long = long / 2
    return (left + middle) / 2


def golden_ratio(sven_, start, direction):
    x0, y0, xd, yd, nd = setup(start, direction, sb=False)
    long = sven_[1] - sven_[0]
    left = sven_[0] + 0.382 * long
    right = sven_[0] + 0.618 * long
    while long > _eps:
        if foo(x0, xd, y0, yd, nd, left) < \
                foo(x0, xd, y0, yd, nd, right):
            sven_[1] = right
            long = sven_[1] - sven_[0]
            left = sven_[0] + 0.382 * long
            right = sven_[0] + 0.618 * long
        elif foo(x0, xd, y0, yd, nd, left) > \
                foo(x0, xd, y0, yd, nd, right):
            sven_[0] = left
            long = sven_[1] - sven_[0]
            left = sven_[0] + 0.382 * long
            right = sven_[0] + 0.618 * long
    return (sven_[0] + sven_[1]) / 2


def dsk(sven_, start, direction):
    x0, y0, xd, yd, nd = setup(start, direction, sb=False)
    dx = abs(sven_[0] - sven_[1])
    f1 = foo(x0, xd, y0, yd, nd, sven_[0])
    f2 = foo(x0, xd, y0, yd, nd, sven_[1])
    f3 = foo(x0, xd, y0, yd, nd, sven_[2])

    x_dsk = sven_[1] + (dx*(f1-f3))/(2*(f1-2*f2+f3))
    f_dsk = foo(x0, xd, y0, yd, nd, x_dsk)
    if (sven_[1]-x_dsk) < _eps and (f2-f_dsk) < _eps:
        return x_dsk
    end = False
    while not end:
        dotes = sorted([sven_[0], sven_[1], sven_[2], x_dsk])
        dsk_index = dotes.index(x_dsk)
        if dsk_index == 0 or dsk_index == 1:
            dotes.remove(dotes[-1])
        elif dsk_index == 2 or dsk_index == 3:
            dotes.remove(dotes[0])
        f1 = foo(x0, xd, y0, yd, nd, dotes[0])
        f2 = foo(x0, xd, y0, yd, nd, dotes[1])
        f3 = foo(x0, xd, y0, yd, nd, dotes[2])
        a1 = (f2-f1)/(dotes[1]-dotes[0])
        a2 = (1/(dotes[2]-dotes[1]))*((f3-f1)/(dotes[2]-dotes[0])-a1)
        x_dsk = (dotes[0]+dotes[1])/2-a1/(2*a2)
        f_dsk = foo(x0, xd, y0, yd, nd, x_dsk)
        tmp_f = [f1, f2, f3, f_dsk]
        tmp_dotes = [dotes[0], dotes[1], dotes[2], x_dsk]
        index_f_min = tmp_f.index(min(tmp_f))
        x_min = tmp_dotes[index_f_min]

        if abs(min(tmp_f)-f_dsk) < _eps and abs(x_min - x_dsk) < _eps:
            end = True
    return x_dsk


def powell(start):
    with open("rgr Zinchuk.txt", "w+") as file:
        file.write(f"Find minimum of function f(x) = 5(x1-8)**2 + x1x2 + 3x2**2\n")
        file.write(f"Minimum reached at X = [8.136, -1.355] and value equal 0\n\n\n")
        i = 0
        position = [start]
        end_ = 0
        while end_ < 4:
            end_ += 1
            file.write(f"--Iteration {i+1}--\n")
            sven_ = sven(position[i], _direction[i])
            file.write(f" Direction = {_direction[i]}, point = {position[i]}\n")
            file.write(f"Lambda interval = {sven_}\n")
            if i == 0:
                lambda_optimal = dichotomy(sven_, position[i], _direction[i])
                file.write(f"Dichotomy method for optimal lambda. Lambda Value = {lambda_optimal}\n")
            elif i == 1:
                lambda_optimal = golden_ratio(sven_, position[i], _direction[i])
                file.write(f"Golden ratio method for optimal lambda. Lambda Value = {lambda_optimal}\n")
            else:
                lambda_optimal = dsk(sven(position[i], _direction[i], dsk_=True), position[i], _direction[i])
                file.write(f"DSK method for optimal lambda. Lambda Value = {lambda_optimal}\n")
            new_x = [0, 0]
            new_x[0] = position[i][0] + lambda_optimal * _direction[i][0]
            new_x[1] = position[i][1] + lambda_optimal * _direction[i][1]
            position.append(new_x)
            file.write(f"New dote {new_x}\n\n\n")

            if (i+1) % 3 == 0:
                old = _direction[i]
                x = position[i+1][0] - position[i-1][0]
                y = position[i+1][1] - position[i-1][1]
                new_dir = [x, y]
                normalized_new = (x**2+y**2) ** 0.5
                new_dir[0] = new_dir[0]/normalized_new
                new_dir[1] = new_dir[1]/normalized_new
                _direction.append(new_dir)
                _direction.append(old)
                _direction.append(new_dir)

            i += 1
            norm_new_x = (position[i][0]**2 + position[i][1]**2) ** 0.5
            prev_norm_x = (position[i-1][0]**2 + position[i-1][1]**2) ** 0.5
            if abs(norm_new_x - prev_norm_x)/norm_new_x <= _eps:
                end_ = 5
        file.write(f"\nEnd of calculation via reaching criteria of ending")
        file.close()
        position[-1][0] = round(position[-1][0], 4)
        position[-1][1] = round(position[-1][1], 4)
        return position


def plot_():
    start, stop, n_values = -5, 15, 2000

    x_values = np.linspace(start, stop, n_values)
    y_values = np.linspace(start, stop, n_values)
    _x, _y = np.meshgrid(x_values, y_values)

    _z = (5 * (_x - 5) ** 2 + _x * _y + 3 * _y ** 2)

    plt.contour(_x, _y, _z)

    path = powell([14.4, 14.4])

    x_coords = []
    for i in range(len(path)):
        x_coords.append(path[i][0])

    y_coords = []
    for i in range(len(path)):
        y_coords.append(path[i][1])

    plt.scatter(x_coords, y_coords)
    plt.plot(x_coords, y_coords)

    plt.show()


if __name__ == '__main__':
    plot_()
