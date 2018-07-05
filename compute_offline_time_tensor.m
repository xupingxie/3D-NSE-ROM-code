cput_start = tic;
[T] = generate_matrices3d_Tenssor(POD_u, POD_v, POD_w, x, e_conn);
cputime = toc(cput_start);
