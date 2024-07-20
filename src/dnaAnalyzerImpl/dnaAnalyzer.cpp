#include "dnaAnalyzer.h"

#include <fstream>
#include <map>
#include <sstream>

#define CONST_D 42 // константа для dnaSearch
#define Q 13       // простое число для dnaSearch

dnaAnalyzer::dnaAnalyzer() = default;

static std::tuple<std::string, std::string>
readFileLines(std::string &filename) {
  std::ifstream file(filename);
  if (!file.good()) {
    file.clear();
    throw std::invalid_argument("File " + filename + " could not be opened.");
  }
  std::string sequence1, sequence2;
  file >> sequence1 >> sequence2;
  file.close();

  return std::make_tuple(sequence1, sequence2);
}

static std::tuple<int, int, int, std::string, std::string>
readFileScoreLines(std::string &filename) {
  std::ifstream file(filename);
  if (!file.good()) {
    file.clear();
    throw std::invalid_argument("Файл " + filename + " не удалось открыть.");
  }

  // Читаем стоимость совпадения, несовпадения и гэпа (если есть)
  int matchScore, mismatchScore, gapScore;
  file >> matchScore >> mismatchScore >> gapScore;

  std::string sequence1, sequence2;
  file >> sequence1 >> sequence2;
  file.close();

  return std::make_tuple(matchScore, mismatchScore, gapScore, sequence1,
                         sequence2);
}

// RabinKarp
std::vector<int> dnaAnalyzer::dnaSearch(std::string &firstFilename,
                                        std::string &secondFilename) {
  std::vector<int> result;
  auto parameters1 = readFileLines(firstFilename);
  auto parameters2 = readFileLines(secondFilename);

  std::string text = std::get<0>(parameters1);
  std::string pattern = std::get<0>(parameters2);

  if (text.empty() || pattern.empty()) {
    throw std::invalid_argument("One of the sequences is empty.");
  }

  int lengthText = text.length();
  int lengthPattern = pattern.length();

  int h = 1;
  for (int i = 0; i < lengthPattern - 1; i++) {
    h = (h * CONST_D) % Q;
  }

  int hashPattern = 0;
  int hashText = 0;

  // считаем хэш для паттерна и первого окна текста
  for (int i = 0; i < lengthPattern; i++) {
    hashText = (CONST_D * hashText + text[i]) % Q;
    hashPattern = (CONST_D * hashPattern + pattern[i]) % Q;
  }

  // перебираем все окна текста
  for (int i = 0; i <= lengthText - lengthPattern; i++) {
    // если хэши совпали, проверяем на совпадение посимвольно
    if (hashPattern == hashText) {
      bool isFound = true;
      for (int j = 0; j < lengthPattern; j++) {
        if (text[i + j] != pattern[j]) {
          isFound = false;
          break;
        }
      }
      if (isFound) {
        result.push_back(i);
      }
    }

    // считаем хэш для следующего окна текста
    if (i < lengthText - lengthPattern) {
      hashText =
          (CONST_D * (hashText - text[i] * h) + text[i + lengthPattern]) % Q;
      if (hashText < 0) {
        hashText = (hashText + Q);
      }
    }
  }

  return result;
}

static void fillMatrix(std::vector<std::vector<long long>> &dp,
                       const std::string &sequence1,
                       const std::string &sequence2, int matchScore,
                       int mismatchScore, int gapScore) {
  int m = sequence1.length();
  int n = sequence2.length();

  // инициализируем первую строку и столбец
  for (int i = 1; i <= m; ++i) {
    dp[i][0] = dp[i - 1][0] + gapScore;
  }

  for (int j = 1; j <= n; ++j) {
    dp[0][j] = dp[0][j - 1] + gapScore;
  }

  // заполняем матрицу на основе рекуррентного соотношения
  for (int i = 1; i <= m; ++i) {
    for (int j = 1; j <= n; ++j) {
      long long match =
          dp[i - 1][j - 1] +
          (sequence1[i - 1] == sequence2[j - 1] ? matchScore : mismatchScore);
      long long deleteScore = dp[i - 1][j] + gapScore;
      long long insertScore = dp[i][j - 1] + gapScore;

      dp[i][j] = std::max({match, deleteScore, insertScore});
    }
  }
}

// NeedlemanWunsh
long long dnaAnalyzer::sequenceAlignmentOptimalScore(std::string &filename) {
  auto parameters = readFileScoreLines(filename);

  int matchScore = std::get<0>(parameters);
  int mismatchScore = std::get<1>(parameters);
  int gapScore = std::get<2>(parameters);
  std::string sequence1 = std::get<3>(parameters);
  std::string sequence2 = std::get<4>(parameters);

  int m = sequence1.length();
  int n = sequence2.length();

  // создаем квадратную матрицу для хранения оценок
  std::vector<std::vector<long long>> dp(m + 1,
                                         std::vector<long long>(n + 1, 0));

  // заполняем матрицу
  fillMatrix(dp, sequence1, sequence2, matchScore, mismatchScore, gapScore);

  // возвращаем оптимальный результат
  return dp[m][n];
}

std::string
dnaAnalyzer::recoveringOptimalSequenceAlignment(std::string &filename) {
  auto parameters = readFileScoreLines(filename);

  int matchScore = std::get<0>(parameters);
  int mismatchScore = std::get<1>(parameters);
  int gapScore = std::get<2>(parameters);
  std::string sequence1 = std::get<3>(parameters);
  std::string sequence2 = std::get<4>(parameters);

  int m = sequence1.length();
  int n = sequence2.length();

  // создаем квадратную матрицу для хранения оценок
  std::vector<std::vector<long long>> dp(m + 1,
                                         std::vector<long long>(n + 1, 0));

  // заполняем матрицу
  fillMatrix(dp, sequence1, sequence2, matchScore, mismatchScore, gapScore);

  // осуществляем обратный проход для восстановления выравнивания
  int i = m, j = n;
  std::string alignedSequence1, alignedSequence2;
  std::string matching;

  while (i > 0 || j > 0) {
    if (i > 0 && j > 0 &&
        dp[i][j] == dp[i - 1][j - 1] + (sequence1[i - 1] == sequence2[j - 1]
                                            ? matchScore
                                            : mismatchScore)) {
      alignedSequence1 = sequence1[i - 1] + alignedSequence1;
      alignedSequence2 = sequence2[j - 1] + alignedSequence2;
      i--;
      j--;
    } else if (i > 0 && dp[i][j] == dp[i - 1][j] + gapScore) {
      alignedSequence1 = sequence1[i - 1] + alignedSequence1;
      alignedSequence2 = '-' + alignedSequence2;
      i--;
    } else {
      alignedSequence1 = '-' + alignedSequence1;
      alignedSequence2 = sequence2[j - 1] + alignedSequence2;
      j--;
    }
    if (alignedSequence1[0] == alignedSequence2[0]) {
      matching += '|';
    } else {
      matching += ' ';
    }
  }
  std::reverse(matching.begin(), matching.end());
  // формируем результат выравнивания с переносами строки
  std::string result;
  result += alignedSequence1 + '\n';
  result += matching + '\n';
  result += alignedSequence2;

  return result;
}

static bool isValidChar(char c) {
  return c == 'A' || c == 'C' || c == 'G' || c == 'T' || c == '.' || c == '?' ||
         c == '+' || c == '*';
}

static std::string processPattern(const std::string &pattern,
                                  std::map<int, char> &plus) {
  std::string processedPattern = "";
  int i = 0;
  while (pattern[i] != '\0') {
    if (isValidChar(pattern[i]) && pattern[i + 1] == '+') {
      processedPattern += '+';
      plus[processedPattern.size() - 1] = pattern[i];
      i++;
    } else if (pattern[i] == '*') {
      processedPattern += '*';
      while (pattern[i + 1] == '*') {
        i++;
      }
    } else if (pattern[i] == '+') {
      throw std::invalid_argument("Invalid pattern");
    } else if (!isValidChar(pattern[i])) {
      throw std::invalid_argument("Invalid character in pattern");
    } else {
      processedPattern += pattern[i];
    }
    i++;
  }
  return processedPattern;
}
void fillMatrix(const std::string &sequence,
                const std::string &processedPattern,
                std::vector<std::vector<bool>> &matrix,
                const std::map<int, char> &plus) {
  matrix[0][0] = true;
  for (int j = 1; j <= processedPattern.size(); j++) {
    if (processedPattern[j - 1] == '*' || processedPattern[j - 1] == '+' ||
        processedPattern[j - 1] == '?') {
      matrix[0][j] = matrix[0][j - 1];
    }
  }

  for (int i = 1; i <= sequence.size(); i++) {
    if (!isValidChar(sequence[i - 1])) {
      throw std::invalid_argument("Invalid character in sequence");
    }

    for (int j = 1; j <= processedPattern.size(); j++) {
      if (processedPattern[j - 1] == '.' ||
          sequence[i - 1] == processedPattern[j - 1]) {
        matrix[i][j] = matrix[i - 1][j - 1];
      } else if (processedPattern[j - 1] == '*') {
        matrix[i][j] = matrix[i - 1][j] || matrix[i][j - 1];
      } else if (processedPattern[j - 1] == '?') {
        matrix[i][j] = matrix[i - 1][j - 1] || matrix[i][j - 1];
      } else if (processedPattern[j - 1] == '+') {
        matrix[i][j] =
            (sequence[i - 1] == plus.at(j - 1)) && matrix[i - 1][j] ||
            matrix[i][j - 1];
      }
    }
  }
}

bool dnaAnalyzer::matchingRegularExpressions(std::string &filename) {
  auto parameters = readFileLines(filename);

  std::string sequence = std::get<0>(parameters);
  std::string pattern = std::get<1>(parameters);

  std::map<int, char> plus;
  std::string processedPattern = processPattern(pattern, plus);

  std::vector<std::vector<bool>> matrix(
      sequence.size() + 1,
      std::vector<bool>(processedPattern.size() + 1, false));
  fillMatrix(sequence, processedPattern, matrix, plus);

  return matrix[sequence.size()][processedPattern.size()];
}

static int kSimilarity(std::string string1, std::string string2) {
  for (int i = 0; i < string1.size(); i++) {
    if (string1[i] == string2[i]) {
      continue;
    }
    std::vector<int> matches;
    for (int j = i + 1; j < string1.size(); j++) {
      if (string1[i] == string2[j] && string1[j] != string2[j]) {
        matches.push_back(j);
        if (string1[j] == string2[i]) {
          std::swap(string2[j], string2[i]);
          return 1 + kSimilarity(string1.substr(i + 1), string2.substr(i + 1));
        }
      }
    }
    int k = string1.size() - 1;
    for (int j : matches) {
      std::swap(string2[i], string2[j]);
      k = std::min(
          k, 1 + kSimilarity(string1.substr(i + 1), string2.substr(i + 1)));
      std::swap(string2[i], string2[j]);
    }
    return k;
  }
  return 0;
}

int dnaAnalyzer::kSimilarityOfTwoSequences(std::string &filename) {
  auto parameters = readFileLines(filename);

  std::string string1 = std::get<0>(parameters);
  std::string string2 = std::get<1>(parameters);

  std::string s1 = string1;
  std::string s2 = string2;
  std::sort(s1.begin(), s1.end());
  std::sort(s2.begin(), s2.end());
  if (s1 != s2) {
    throw std::invalid_argument("Strings " + string1 + " and " + string2 +
                                " are not anagrams.");
  }

  return kSimilarity(string1, string2);
}

std::string dnaAnalyzer::minimumWindowSubstring(std::string &filename) {
  auto parameters = readFileLines(filename);

  std::string text = std::get<0>(parameters);
  std::string pattern = std::get<1>(parameters);

  std::map<char, int> alphabet = {{'A', 0}, {'G', 0}, {'C', 0}, {'T', 0}};
  for (auto &c : pattern) {
    alphabet[c]++;
  }

  int count = pattern.size();
  int minLen = INT_MAX;
  int head = 0;
  int currHead = 0;

  for (int i = 0; i < text.size(); i++) {
    if (--alphabet[text[i]] >= 0) {
      count--;
    }
    while (count == 0) {
      if (minLen > i - currHead + 1) {
        minLen = i - currHead + 1;
        head = currHead;
      }
      if (++alphabet[text[currHead]] > 0) {
        count++;
      }
      currHead++;
    }
  }

  if (minLen == INT_MAX) {
    return "";
  } else {
    return text.substr(head, minLen);
  }
}

dnaAnalyzer::~dnaAnalyzer() = default;
