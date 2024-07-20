#include "ConsoleApp.h"

#include <sys/resource.h>

#include <iostream>
#include <limits>

#include "StringUtils.h"

void ConsoleApp::askToPressAnyKey() {
  std::cout << "Press ENTER to continue . . . ";
  std::cin.ignore();
  std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
}

bool ConsoleApp::tryParseIntFromInput(int *pInt) {
  std::string input;
  std::cin >> input;
  StringUtils::trim(input);

  try {
    *pInt = std::stoi(input);
  } catch (const std::exception &e) {
    return false;
  }

  return true;
}

command ConsoleApp::getCommand() {
  while (true) {
    std::cout
        << std::endl
        << "\n_____________________________________________________________\n"
           "Please, choose an action, then hit a number and press 'enter': \n"
           "#1: exact DNA search (Rabin-Karp algorithm)\n"
           "#2: sequence alignment: calculation of the optimal score "
           "(Needleman-Wunsch algorithm)\n"
           "#3: recovering optimal sequence alignment (Needleman-Wunsch "
           "algorithm)\n"
           "#4: check if a sequence matches a regular expression\n"
           "#5: calculation the k-similarity of two sequences\n"
           "#6: minimum window substring for a sequence (BONUS)\n"
           "#7: quit the program\n"
           "_____________________________________________________________"
        << std::endl;

    int number;
    if (!tryParseIntFromInput(&number)) {
      std::cout << "Invalid input. Input should be of integer type. Try again."
                << std::endl;
      askToPressAnyKey();
      continue;
    }

    if (number < MIN_COMMAND || number > MAX_COMMAND) {
      std::cout << number
                << " is out of bounds of expected command numbers. Try again"
                << std::endl;
      askToPressAnyKey();
      continue;
    }

    return (enum command)number;
  }
}

void ConsoleApp::executeCommand(command command) {
  dnaAnalyzer dnaAnalyzer;
  switch (command) {
  case dnaSearch: {
    std::cout << "Please, provide two filenames containing the sequences: "
              << std::endl;

    std::string firstFilename;
    std::string secondFilename;
    std::cin >> firstFilename >> secondFilename;

    auto start =
        std::chrono::high_resolution_clock::now(); // Замеряем начальное время
    struct rusage startUsage;
    getrusage(RUSAGE_SELF, &startUsage);
    std::vector<int> result =
        dnaAnalyzer::dnaSearch(firstFilename, secondFilename);
    struct rusage endUsage;
    getrusage(RUSAGE_SELF, &endUsage);
    auto end =
        std::chrono::high_resolution_clock::now(); // Замеряем конечное время
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(
        end - start); // Получаем разницу в микросекундах
    long maxResidentSetSize = endUsage.ru_maxrss / 1024;

    std::cout << "result:\n";
    if (result.empty()) {
      std::cout << "no matches found";
    } else {
      for (int i : result) {
        std::cout << i << " ";
      }
    }
    std::cout << std::endl;
    std::cout << "Elapsed (wall clock) time: " << duration.count()
              << " milliseconds" << std::endl;
    std::cout << "Max Resident Set Size: " << maxResidentSetSize
              << " kilobytes\n";
    break;
  }
  case sequenceAlignmentOptimalScore: {
    std::cout << "Please, provide the filename containing\n"
                 "match score, mismatch, and gap;\n"
                 "two sequences to align: "
              << std::endl;
    std::string filename;
    std::cin >> filename;
    auto start =
        std::chrono::high_resolution_clock::now(); // Замеряем начальное время
    struct rusage startUsage;
    getrusage(RUSAGE_SELF, &startUsage);
    long long result = dnaAnalyzer::sequenceAlignmentOptimalScore(filename);
    struct rusage endUsage;
    getrusage(RUSAGE_SELF, &endUsage);
    auto end =
        std::chrono::high_resolution_clock::now(); // Замеряем конечное время
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(
        end - start); // Получаем разницу в микросекундах
    long maxResidentSetSize = endUsage.ru_maxrss / 1024;

    std::cout << "result:\n" << result << std::endl;
    std::cout << "Elapsed (wall clock) time: " << duration.count()
              << " milliseconds" << std::endl;
    std::cout << "Max Resident Set Size: " << maxResidentSetSize
              << " kilobytes\n";
    break;
  }
  case recoveringOptimalSequenceAlignment: {
    std::cout << "Please, provide the filename containing\n"
                 "match score, mismatch, and gap;\n"
                 "two sequences to align: "
              << std::endl;
    std::string filename;
    std::cin >> filename;
    auto start =
        std::chrono::high_resolution_clock::now(); // Замеряем начальное время
    struct rusage startUsage;
    getrusage(RUSAGE_SELF, &startUsage);
    std::string result =
        dnaAnalyzer::recoveringOptimalSequenceAlignment(filename);
    struct rusage endUsage;
    getrusage(RUSAGE_SELF, &endUsage);
    auto end =
        std::chrono::high_resolution_clock::now(); // Замеряем конечное время
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(
        end - start); // Получаем разницу в микросекундах
    long maxResidentSetSize = endUsage.ru_maxrss / 1024;

    std::cout << "result:" << std::endl;
    std::cout << dnaAnalyzer::sequenceAlignmentOptimalScore(filename)
              << std::endl;
    std::cout << result << std::endl;
    std::cout << "Elapsed (wall clock) time: " << duration.count()
              << " milliseconds" << std::endl;
    std::cout << "Max Resident Set Size: " << maxResidentSetSize
              << " kilobytes\n";
    break;
  }
  case matchingRegularExpressions: {
    std::cout << "Please, provide the filename containing\n"
                 "a sequence and a regular expression: "
              << std::endl;
    std::string filename;
    std::cin >> filename;
    auto start =
        std::chrono::high_resolution_clock::now(); // Замеряем начальное время
    struct rusage startUsage;
    getrusage(RUSAGE_SELF, &startUsage);
    bool result = dnaAnalyzer::matchingRegularExpressions(filename);
    struct rusage endUsage;
    getrusage(RUSAGE_SELF, &endUsage);
    auto end =
        std::chrono::high_resolution_clock::now(); // Замеряем конечное время
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(
        end - start); // Получаем разницу в микросекундах
    long maxResidentSetSize = endUsage.ru_maxrss / 1024;

    std::cout << "result:\n" << (result ? "True" : "False") << std::endl;
    std::cout << "Elapsed (wall clock) time: " << duration.count()
              << " milliseconds" << std::endl;
    std::cout << "Max Resident Set Size: " << maxResidentSetSize
              << " kilobytes\n";
    break;
  }
  case kSimilarityOfTwoSequences: {
    std::cout << "Please, provide the filename containing two sequences: "
              << std::endl;
    std::string filename;
    std::cin >> filename;
    auto start =
        std::chrono::high_resolution_clock::now(); // Замеряем начальное время
    struct rusage startUsage;
    getrusage(RUSAGE_SELF, &startUsage);
    int result = dnaAnalyzer::kSimilarityOfTwoSequences(filename);
    struct rusage endUsage;
    getrusage(RUSAGE_SELF, &endUsage);
    auto end =
        std::chrono::high_resolution_clock::now(); // Замеряем конечное время
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(
        end - start); // Получаем разницу в микросекундах
    long maxResidentSetSize = endUsage.ru_maxrss / 1024;

    std::cout << "result:\n" << result << std::endl;
    std::cout << "Elapsed (wall clock) time: " << duration.count()
              << " milliseconds" << std::endl;
    std::cout << "Max Resident Set Size: " << maxResidentSetSize
              << " kilobytes\n";

    break;
  }
  case minimumWindowSubstring: {
    std::cout << "Please, provide the filename containing\n"
                 "a sequence and a substring: "
              << std::endl;
    std::string filename;
    std::cin >> filename;
    auto start =
        std::chrono::high_resolution_clock::now(); // Замеряем начальное время
    struct rusage startUsage;
    getrusage(RUSAGE_SELF, &startUsage);
    std::string result = dnaAnalyzer::minimumWindowSubstring(filename);
    struct rusage endUsage;
    getrusage(RUSAGE_SELF, &endUsage);
    auto end =
        std::chrono::high_resolution_clock::now(); // Замеряем конечное время
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(
        end - start); // Получаем разницу в микросекундах
    long maxResidentSetSize = endUsage.ru_maxrss / 1024;

    std::cout << "result:\n" << result << std::endl;
    std::cout << "Elapsed (wall clock) time: " << duration.count()
              << " milliseconds" << std::endl;
    std::cout << "Max Resident Set Size: " << maxResidentSetSize
              << " kilobytes\n";
    break;
  }
  case quit: {
    break;
  }
  }
}

void ConsoleApp::start() {

  while (true) {
    command command = getCommand();
    if (command == quit) {
      return;
    }
    try {
      executeCommand(command);
      askToPressAnyKey();
    } catch (std::invalid_argument &ex) {
      std::cout << "invalid argument: " << ex.what() << std::endl;
      askToPressAnyKey();
    } catch (std::exception &ex) {
      std::cout << "unexpected exception: " << ex.what() << std::endl;
      askToPressAnyKey();
    }
  }
}
