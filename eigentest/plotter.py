import sys
import numpy as np
import matplotlib.pyplot as plt

if __name__ == '__main__':
    filename = sys.argv[1]

    lines = []
    nres = None
    with open(filename) as fp:
        for l in fp:
            ls = l.strip()
            if len(ls) == 0 or ls[0] == '#':
                continue

            ls = ls.split()

            ls[0] = int(ls[0])
            ls[1:] = map(float, ls[1:])
            lines.append(ls)


    nres = (len(lines[0]) - 1) // 4
    if nres == 3:
        lapack_load, lapack_trans, lapack_comp, lapack_error = (None,)*4
        (
            dim,
            cusolver_load , cusolver_trans , cusolver_comp , cusolver_error ,
            magma_cpu_load, magma_cpu_trans, magma_cpu_comp, magma_cpu_error,
            magma_gpu_load, magma_gpu_trans, magma_gpu_comp, magma_gpu_error
        ) = zip(*lines)
    elif nres == 4:
        (
            dim,
            lapack_load   , lapack_trans   , lapack_comp   , lapack_error   ,
            cusolver_load , cusolver_trans , cusolver_comp , cusolver_error ,
            magma_cpu_load, magma_cpu_trans, magma_cpu_comp, magma_cpu_error,
            magma_gpu_load, magma_gpu_trans, magma_gpu_comp, magma_gpu_error
        ) = zip(*lines)

    dim = np.array(dim)
    for label, color, load, trans, comp, error in [('Lapack'   , 'k', lapack_load   , lapack_trans   , lapack_comp   , lapack_error   ),
                                                   ('CuSolver' , 'r', cusolver_load , cusolver_trans , cusolver_comp , cusolver_error ),
                                                   ('Magma/cpu', 'g', magma_cpu_load, magma_cpu_trans, magma_cpu_comp, magma_cpu_error),
                                                   ('Magma/gpu', 'b', magma_gpu_load, magma_gpu_trans, magma_gpu_comp, magma_gpu_error)]:
        if load is None: continue
        load, trans, comp, error = tuple(np.array(x) for x in (load, trans, comp, error))
        filt = error < 1.0e-5
        plt.plot(dim[filt], comp[filt]                       , label=f'{label} (comp)' , ls='--' , color=color, marker='+')
        plt.plot(dim[filt], load[filt]                       , label=f'{label} (load)' , ls='-.' , color=color, marker='x')
        #plt.plot(dim[filt], trans[filt]                      , label=f'{label} (trans)', ls=':', color=color, marker='*')
        plt.plot(dim[filt], trans[filt]+comp[filt]           , label=f'{label} (comp+trans)', ls='-'  , color=color, marker='.')
        #plt.plot(dim[filt], load[filt]+trans[filt]+comp[filt], label=f'{label} (total)', ls='-'  , color=color, marker='.')

    plt.xlim(left=0)
    plt.ylim(bottom=0)
    plt.xlabel('Matrix dimension')
    plt.ylabel('Realtime [sec]')
    plt.legend()
    plt.show()

    if lapack_comp is None:
        sys.exit(0)

    # Speedup
    base_comp = np.array(lapack_comp)
    for label, color, trans, comp, error in [('Lapack'   , 'k', lapack_trans   , lapack_comp   , lapack_error   ),
                                             ('CuSolver' , 'r', cusolver_trans , cusolver_comp , cusolver_error ),
                                             ('Magma/cpu', 'g', magma_cpu_trans, magma_cpu_comp, magma_cpu_error),
                                             ('Magma/gpu', 'b', magma_gpu_trans, magma_gpu_comp, magma_gpu_error)]:
        if load is None: continue

        trans, comp, error = tuple(np.array(x) for x in (trans, comp, error))

        filt = error < 1.0e-5
        plt.plot(dim[filt], base_comp[filt]/comp[filt]              , label=f'{label} (comp)' , ls='--' , color=color, marker='+')
        plt.plot(dim[filt], base_comp[filt]/(trans[filt]+comp[filt]), label=f'{label} (comp+trans)', ls='-'  , color=color, marker='.')


    plt.xlim(left=0)
    plt.ylim(bottom=0)
    plt.xlabel('Matrix dimension')
    plt.ylabel('Speedup')
    plt.legend()
    plt.show()
