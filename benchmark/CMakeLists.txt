set(benchmarks
  smartcomponents
  )

foreach(benchmark ${benchmarks})
  add_executable(${benchmark} ${benchmark}.cpp)
  target_link_libraries(${benchmark} ${Helium_LIBRARIES})
endforeach(benchmark ${benchmarks})
