#!/bin/python3
import numpy as np

def grid_convergence_index(
        # area_or_volume,
        num_cells: tuple[float, float, float],
        phi: tuple[float, float, float],
        dimension='2d',
    ):

    # description
    """
    I. B. Celik, U. Ghia, P. J. Roache, C. J. Freitas, H. Coleman, und P. E. Raad,
    „Procedure for Estimation and Reporting of Uncertainty
        Due to Discretization in CFD Applications“,
    Journal of Fluids Engineering, Bd. 130, Nr. 078001, Juli 2008, doi: 10.1115/1.2960953.
    each parameter needs to be ordered from the attributes of the finest mesh to the coarsest mesh
    # :param area_or_volume: (1,) number array
        representing the whole area (2d) or volume (3d) of the domain
    :type num_cells:
        * (3,) number array representing the cells count for each mesh, orderd in
        * sorted from finest to coarsest mesh / biggest cell count to least cell count
    :param phi: (3,) number array representing the solution for each mesh
    :type dimension: '2d' or '3d'
    """

    # step 1
    print("It is assumed that the area, respectively the volume, remained the same "
          "in all three studies!")
    match dimension:
        case '2d':
            # h: (3,) = np.sqrt(area_or_volume / num_cells)
            r21 = np.sqrt( np.divide(num_cells[0], num_cells[1]) )
            r32 = np.sqrt( np.divide(num_cells[1], num_cells[2]) )
        case '3d':
            # h: (3,) = np.cbrt(area_or_volume / num_cells)
            r21 = np.cbrt( np.divide(num_cells[0], num_cells[1]) )
            r32 = np.cbrt( np.divide(num_cells[1], num_cells[2]) )
        case _:
            print("Key 'dimension' only has Values: '2d', '3d'.")
            exit(1)
    # r21 = h[1] / h[0]
    # r32 = h[2] / h[1]

    # step 3
    if (num_cells[2] > num_cells[1]) or (num_cells[1] > num_cells[0]):
        print("It seems like the number of cells are not in correct order!\n"
              "The attributes of the studies need to be sorted from finest mesh to coarsest mesh.\n"
              "The results are still going to be calculated.")
    epsilon21 = np.subtract(phi[1], phi[0])
    epsilon32 = np.subtract(phi[2], phi[1])
    if (np.round(epsilon21, 3) == 0) or (np.round(epsilon32, 3) == 0):
        print(f"The value for epsilon21 is {epsilon21} "
              f"which is approximately {np.round(epsilon32, 3)}.\n"
              f"The value for epsilon32 is {epsilon32} "
              f"which is approximately {np.round(epsilon32, 3)}.\n"
              "One of the values is too close to zero, which is why the formulas applied here "
              "do not work.")
        exit(1)
    q = 0  # initial guess
    p = np.divide(
            np.abs(np.add(
                    np.log(np.abs(
                        np.divide(epsilon32, epsilon21)
                    )),
                    q
                )),
            np.log(r21)
        )
    precison_r = 2
    if np.round(r21, precison_r) != np.round(r32, precison_r):
        s = 1 * np.sign( np.divide(epsilon32, epsilon21) )
        p_iterative = [0]
        precision_p = 5
        print("A refinement factor between the studies did not remained constant.\n"
              f"A fixed point iteration is used to calculate the final p up to "
              f"{precision_p} decimal places.")
        while np.round(p, precision_p) != np.round(p_iterative[-1], precision_p):
            print(f"q = {q}")
            print(f"p = {p}")
            p_iterative.append(p)
            q = np.log(np.divide(
                        np.subtract(np.power(r21, p), s),
                        np.subtract(np.power(r32, p), s)
                ))
            p = np.divide(
                    np.abs(np.add(
                            np.log(np.abs(
                                np.divide(epsilon32, epsilon21)
                            )),
                            q
                        )),
                    np.log(r21)
                )
        del precision_p

    # step 4
    phi21_extrapolated = \
        np.divide(
            np.subtract(
                np.multiply(np.power(r21, p), phi[0]),
                phi[1]
            ),
            np.subtract(np.power(r21, p), 1)
        )
    phi32_extrapolated = \
        np.divide(
            np.subtract(
                np.multiply(np.power(r32, p), phi[1]),
                phi[2]
            ),
            np.subtract(np.power(r32, p), 1)
        )

    # step 5
    error21_approximate = \
        np.abs(np.divide(
                np.subtract(phi[0], phi[1]),
                phi[0]
            ))
    error21_extrapolated = \
        np.abs(np.divide(
                np.subtract(phi21_extrapolated, phi[0]),
                phi21_extrapolated
            ))
    gci21_fine = \
        np.divide(
            np.multiply(1.25, error21_approximate),
            np.subtract(np.power(r21, p), 1)
        )

    # print(f"h = {h}")
    print(f"r21 = {r21}")
    print(f"r32 = {r32}")
    print(f"epsilon21 = {epsilon21}")
    print(f"epsilon32 = {epsilon32}")
    print(f"p = {p}")
    print(f"phi21_extrapolated = {phi21_extrapolated}")
    print(f"phi32_extrapolated = {phi32_extrapolated}")
    print(f"error21_approximate = {error21_approximate}")
    print(f"error21_extrapolated = {error21_extrapolated}")
    print(f"gci21_fine = {gci21_fine}")

    return
