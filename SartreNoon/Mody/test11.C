#include <iostream>
#include <string>

void MyAnalysisFunction() {
    //std::string userInput;
    double a ;
    // Ask for input
    std::cout << "Please enter your input: ";
    //std::getline(std::cin, userInput);
    std::cin >> a;

    // Process the input
    std::cout << "You entered: " << a << std::endl;
}

void test11() {
    // Call your analysis function
    MyAnalysisFunction();
}
