import argparse
import numpy as np
from util import *

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-d', '--dir', type=str, help='Path to dir')
    parser.add_argument('-t', '--tau', type=float, default=1e-3,
                        help='Regularization parameter')
    args = parser.parse_args()

    dir = args.dir
    tau = args.tau

    X = np.load(f'{dir}/X.npy')
    Y = np.load(f'{dir}/Y.npy')
    loading_init = effectiveInit(X, Y)
    loading, _ = twoDcca(X, Y, loading_init=loading_init, x_regular=tau, y_regular=tau, iter_max=1000)

    np.save(f'{dir}/Lx.npy', loading[0])
    np.save(f'{dir}/Rx.npy', loading[1])
    np.save(f'{dir}/Ly.npy', loading[2])
    np.save(f'{dir}/Ry.npy', loading[3])
