add_test_exec(test_har_ff test_har_ff.cpp)
add_test(NAME test_har_ff COMMAND test_har_ff)

add_test_exec(orca_setup orca_setup.cpp)
add_test(NAME orca_setup COMMAND orca_setup)

set_tests_properties(test_har_ff PROPERTIES FIXTURES_REQUIRED orca)
set_tests_properties(orca_setup  PROPERTIES FIXTURES_SETUP    orca)