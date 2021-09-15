program demo_fftpack_hilbert

    use fftpack, only: hilbert, dp

    print *, hilbert([complex(kind=dp) :: 1, 2, 3, 4   ])
    print *, hilbert([complex(kind=dp) :: 1, 2, 3, 4, 5])

    !! [(1.000,1.000), (2.000,-1.000),  (3.000,-1.000),  (4.000,1.000)]
    !! [(1.000,1.701), (2.000,-1.376), (3.000,-0.6498), (4.000,-1.376), (5.000,1.701)]

end program demo_fftpack_hilbert