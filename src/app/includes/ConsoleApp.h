#ifndef CONSOLE_APP_H
#define CONSOLE_APP_H

#include <string>

#include "dnaAnalyzer.h"

enum command {
  dnaSearch = 1,
  sequenceAlignmentOptimalScore = 2,
  recoveringOptimalSequenceAlignment = 3,
  matchingRegularExpressions = 4,
  kSimilarityOfTwoSequences = 5,
  minimumWindowSubstring = 6,
  quit = 7
};

#define MIN_COMMAND dnaSearch
#define MAX_COMMAND quit

class ConsoleApp {
private:
  static void askToPressAnyKey();
  static bool tryParseIntFromInput(int *pInt);
  static command getCommand();
  static void executeCommand(command command);

public:
  static void start();
};

#endif
