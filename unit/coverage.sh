LVM_PROFILE_FILE="dgt-unit-tests.profraw" ./dgt-unit-tests
/Library/Developer/CommandLineTools/usr/bin/llvm-profdata merge -sparse dgt-unit-tests.profraw -o dgt-unit-tests.profdata
/Library/Developer/CommandLineTools/usr/bin/llvm-cov show ./dgt-unit-tests -instr-profile=dgt-unit-tests.profdata
/Library/Developer/CommandLineTools/usr/bin/llvm-cov report ./dgt-unit-tests -instr-profile=dgt-unit-tests.profdata
