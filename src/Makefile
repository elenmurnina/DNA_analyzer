APP_NAME				= DNAAnalyzer

DNAANALYZER_DIR			= dnaAnalyzerImpl
DNAANALYZER_HEADERS_DIR	= $(DNAANALYZER_DIR)/includes
DNAANALYZER_INCLUDES	= -I$(DNAANALYZER_HEADERS_DIR)
DNAANALYZER_SRC			= $(shell find $(DNAANALYZER_DIR) -name "*.cpp")
DNAANALYZER_OBJ			= $(patsubst %.cpp, %.o, $(DNAANALYZER_SRC))

APP_DIR					= app
APP_HEADERS_DIR			= $(APP_DIR)/includes
APP_INC					= -I$(APP_HEADERS_DIR) $(DNAANALYZER_INCLUDES)
APP_SRC					= $(shell find $(APP_DIR) -name "*.cpp")
APP_OBJ					= $(patsubst %.cpp, %.o, $(APP_SRC))

DATASETS_PATH			= tests/datasets/

TEST					= test
TESTS_DIR				= tests
TEST_SRC				= $(shell find $(TESTS_DIR) -name "*.cpp")

COV_SRC					= $(DNAANALYZER_SRC)

TEST_OBJ				= $(patsubst %.cpp, %.o, $(TEST_SRC))
COV_OBJ					= $(patsubst %.cpp, %.o, $(COV_SRC))
REPORT					= report

CXX						= g++
FLAGS					= -std=c++17 -DDATASETS_PATH=\"$(DATASETS_PATH)\"
TEST_FLAGS				= -lgtest -lgtest_main
GCOV_FLAGS				= --coverage
REPORT_FLAG				= --rc derive_function_end_line=0

app/%.o: app/%.cpp Makefile
	$(CXX) $(FLAGS) -c $(APP_INC) $< -o $@

dnaAnalyzerImpl/%.o: dnaAnalyzerImpl/%.cpp Makefile
	$(CXX) $(FLAGS) -c $(DNAANALYZER_INCLUDES) $< -o $@

tests/%.o: tests/%.cpp Makefile
	$(CXX) $(FLAGS) -c $(APP_INC) $(DNAANALYZER_INCLUDES) $< -o $@

.cpp.o: Makefile
	$(CXX) $(FLAGS) -c $(APP_INC) $(DNAANALYZER_INCLUDES) $< -o $(<:.cpp=.o)

all:			app

app:			$(APP_OBJ) $(DNAANALYZER_OBJ) Makefile
				$(CXX) $(FLAGS) $(APP_OBJ) $(DNAANALYZER_OBJ) -o $(APP_NAME)

tests:			$(TEST_OBJ) $(COV_OBJ) Makefile
				$(CXX) $(FLAGS) $(TEST_OBJ) $(COV_OBJ) $(TEST_FLAGS) -o $(TEST)
				./$(TEST)

clean:
				@rm -f $(DNAANALYZER_OBJ) $(APP_OBJ) $(TEST_OBJ) $(COV_OBJ) $(TEST) \
				**/*.gcno **/*.gcda **/*.info
				@rm -rf ./$(REPORT)

fclean:			clean
				@rm -rf $(APP_NAME)

format:
				cp ../materials/linters/clang-format.txt ./
				clang-format -i $(shell find . -name "*.cpp" -o -name "*.h")
				rm clang-format.txt

style:
				@cp ../materials/linters/clang-format.txt ./
				@clang-format -n $(shell find . -name "*.cpp" -o -name "*.h")
				@rm clang-format.txt

gcov_report:    FLAGS += -fprofile-arcs -ftest-coverage
gcov_report:	$(TEST_OBJ) $(COV_OBJ)
				$(CXX) $(FLAGS) $(TEST_OBJ) $(COV_OBJ) $(TEST_FLAGS) --coverage -o $(TEST)
				./$(TEST)
				lcov -t "$(TEST)" -o gcov.info $(REPORT_FLAG) --no-external -c -d .
				genhtml -o $(REPORT) gcov.info $(REPORT_FLAG)
				@rm -f **/*.gcno **/*.gcda **/*.info

open_report:
				open ./$(REPORT)/index.html

.PHONY:			all app clean fclean tests gcov_report open_report format style