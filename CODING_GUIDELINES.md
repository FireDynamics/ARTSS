# Contributing Guidelines

The coding guidelines of ARTSS are mainly based on the [google coding guidelines](https://google.github.io/styleguide/cppguide.html). Accordingly, our rules are not one-to-one with the google standard, therefore we ask you to have a look at the following rules if you want to contribute to ARTSS. It is not necessary to know the google coding guidelines, it is sufficient to read this file.

## Naming Convention
#### Variable Naming:
  We use underscores. The main reason for this is that we have variables which have different meanings depending on upper and lower case (see class Domain.cpp).

  Example:
  
  `int exemplary_variable_name;`

- Naming Branches: Format: `Issue_XY_short_description`

  Example:
  
  `branch: issue_15_coding_guidelines`
 
  
#### Pull Request:
  - Naming Format: `Issue XY short description`

    Example:
  
    `Issue 15 Coding Guidelines`

  - Naming final commit of PR (to master branch): Format: `Issue XY short description (#PR)`. This will be done automatically by GitHub if the PR is named accordingly.

    Example:
  
    `final commit: Issue 15 Coding Guidelines (#113) [description of changes...]`
   
  - Linking to issue. Either link PR to the issue manually or use the phrase "should close #XY" int the description for automatic linking. 
  - Conclude changes to the PR with a comment summarising the changes to the PR before requesting a review. This is to help the reviewer understand what changes should be part of the PR and also to identify code that may have been added incorrectly. Imagine that you want to avoid having to guess what the other person was trying to achieve with these changes when reviewing someone else's code. This should also be the text that is added as a merge commit message.


## Formatting
#### Parenthesis:
  - the opening parenthesis is in the same line as the opening statement and NOT on a separate line
  - do not drop parenthesis even if the the body of a statement consist of a single line

    Example:
    ```
    if (check) {
       body
    }
    ```
#### Pointers and References:
  We are using the style where you emphasise the type of the pointed-to data. `someType *somePtr` rather than `someType* somePtr`. The reasoning is mainly that the pointer belongs to the data type because if you define multiple pointers in one line, you have to write: `someType *p1, *p2` while `someType* p1, p2` creates only one pointer.

    Example:

    `int *pointer1, *pointer2;`
    
#### Line Break/Hard Wrap
  100 characters per line is the guideline here but not a rule which has to be followed everywhere. We are adding line breaks if it makes the code more readable, which in most cases means were a breaking a line out of semantic reasons. In which case the line can be much shorter than 100 characters. There are also some special cases where it is not possible to break a line.
   - Function names/calls are longer than 100 characters
     For readability the whole argument list will be wrapped to the next line. If the argument lists fits in the next line, it should look like this:
     ```
     void Classname::some_long_function_name(
        int argument1, int argument2) {
     [...]
     }
     ```
     If the argument list is too long to fit on one line, each argument gets its own line or is semantically wrapped. 
     ```
     void Classname::some_long_function_name(
        int argument_long1,
        int argument_long2)
     [...]
     }
     ```
     Sometimes it makes more sense to group the argument list, like this:
     ```
     void Classname::some_long_function_name(
         int argument_long_start1, int argument_long_end1,
         int argument_long_start2, int argument_long_end2)
     [...]
     }
     ```
     Our top priority is to the readability while trying to be as consistent as possible.

## Includes
#### Include What You Use
  If a source or header file refers to a symbol defined elsewhere, the file should directly include a header file which properly intends to provide a declaration or definition of that symbol. It should not include header files for any other reason. 

#### Names and Order of Includes
  Include headers in the following order: Related header, C system headers, C++ standard library headers, other libraries' headers, your project's headers.

  All of a project's header files should be listed as descendants of the project's source directory without use of UNIX directory aliases . (the current directory) or .. (the parent directory). For example, google-awesome-project/src/base/logging.h should be included as:

  #include "base/logging.h"

  In dir/foo.cc or dir/foo_test.cc, whose main purpose is to implement or test the stuff in dir2/foo2.h, order your includes as follows:
  ```
      1. dir2/foo2.h.
      2. A blank line
      3. C system headers (more precisely: headers in angle brackets with the .h extension), e.g., <unistd.h>, <stdlib.h>.
      4. A blank line
      5. C++ standard library headers (without file extension), e.g., <algorithm>, <cstddef>.
      6. A blank line
      7. Other libraries' .h files.
      8. Your project's .h files.
  ```
  Separate each non-empty group with one blank line.

## Header Files
#### Self contained:
  Headers should be self contained, which means it should contain all necessary imports without relying on import of other (imported) headers.

  Example: `Solution.h` uses `std::vector` accordingly it has to import `std::vector`

#### define guards:
  Our include guards have the format ARTSS_DIRNAME_FILENAME_H_

  Example: `Soluction.h` in directory `analysis`
  
  ```
  #ifndef ARTSS_ANALYSIS_SOLUTION_H_
  #define ARTSS_ANALYSIS_SOLUTION_H_
  
  [...]
  
  #endif /* ARTSS_ANALYSIS_SOLUTION_H_ */
  ```

 

## Scoping
#### The auto keyword
  We like to use the `auto` type specifier to improve code readability. For example in situations where the variable type matches the initialiser expression. For more detailed information see Clang Tidy rule `moderinze-use-auto`.
  
#### Static and global variables and singletons
  ARTSS currently uses some static and global variables as well two singletons, please refrain from adding new ones except for a exceptional good reason.

As we are only a small group of developers, many parts of the code are not yet compliant with our code guidelines and will be improved accordingly as we make changes. Accordingly, more rules may be added to our coding guidelines if we determine that it is necessary.

Please feel free to contact us if you have any questions.
