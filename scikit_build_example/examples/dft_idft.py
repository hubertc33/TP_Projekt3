import scikit_build_example as s

a="sin"
b=s.dft(s.gen_signal(a,10))
s.show_plot(b)
s.show_plot(s.idft(b))