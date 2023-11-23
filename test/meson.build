testdrive_dep = dependency('test-drive', fallback: ['test-drive', 'testdrive_dep'])

tests = [
  'fft',
  'rfft',
  'dct',
  'utils'
]

test_srcs = files(
  'test_fftpack.f90'
)
foreach t : tests
  test_srcs += files('test_fftpack_@0@.f90'.format(t.underscorify()))
endforeach

tester = executable(
  'tester',
  sources: test_srcs,
  dependencies: [fftpack_dep, testdrive_dep]
)

foreach t : tests
  test(t, tester, args: t)
endforeach
