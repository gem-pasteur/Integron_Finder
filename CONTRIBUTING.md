# How to contribute

Welcome! By visiting this page,
you've taken the first step toward becoming a contributor to integron_finder!

## What skills do I need?

If you know Python3, you can contribute to the code.
We also use [pandas](https://pandas.pydata.org/) and
[numpy](http://www.numpy.org/) libraries.

And there are also many ways to contribute to the project without programming.
If you'd like to get involved in design,
support, testing, documentation, or other types of contributions.


## How Can I Contribute?

We use Github to host the project, to submit an issue or a pull request,
you have to create a GitHub account.

### Reporting Bugs

This section guides you through submitting a bug report for integron_finder.
Following these guidelines helps maintainers and the community understand your report,
reproduce the behavior, and find related reports.

Before creating bug reports, please check [this list](#before-submitting-a-bug-report)
as you might find out that you don't need to create one.
When you are creating a bug report, please
[include as many details as possible](#how-do-i-submit-a-good-bug-report).
Fill out [the required template](.github/ISSUE_TEMPLATE.md), the information it asks for helps us resolve issues faster.

> **Note:**
> If you find a **Closed** issue that seems like it is the same thing that you're experiencing,
> open a new issue and include a link to the original issue in the body of your new one.

#### Before Submitting A Bug Report

* **Perform a [cursory search](https://github.com/gem-pasteur/Integron_Finder/issues?q=is%3Aopen+is%3Aissue+label%3Abug)** to see if the problem has already been reported.
If it has **and the issue is still open**, add a comment to the existing issue instead of opening a new one.


#### How Do I Submit A (Good) Bug Report?

Bugs are tracked as [GitHub issues](https://guides.github.com/features/issues/).
Create an issue on that repository and provide the following information by filling in the template.
Explain the problem and include additional details to help maintainers reproduce the problem:

* **Use a clear and descriptive title** for the issue to identify the problem.
* **Describe the exact steps which reproduce the problem** in as many details as possible.
  Which command exactly you used in the terminal.
  When listing steps, **don't just say what you did, but explain how you did it**.
* **Provide specific examples to demonstrate the steps**.
  Include links to files or GitHub projects, or copy/pasteable snippets, which you use in those examples.
  If you're providing snippets in the issue, use [Markdown code blocks](https://help.github.com/articles/markdown-basics/#multiple-lines).
* **Describe the behavior you observed after following the steps** and point out what exactly is the problem with that behavior.
* **Explain which behavior you expected to see instead and why.**
* **If you're reporting that integron_finder crashed**,
  include a crash report with a stack trace from the operating system.
  Include the crash report in the issue in a [code block](https://help.github.com/articles/markdown-basics/#multiple-lines),
  a [file attachment](https://help.github.com/articles/file-attachments-on-issues-and-pull-requests/),
  or put it in a [gist](https://gist.github.com/) and provide link to that gist.
* **If the problem wasn't triggered by a specific action**, describe what you were doing before the problem happened
  and share more information using the guidelines below.

##### Provide more context by answering these questions:

* **Did the problem start happening recently** (e.g. after updating to a new version of integron_finder) or was this always a problem?
* If the problem started happening recently, **can you reproduce the problem in an older version of integron_finder?**
  What's the most recent version in which the problem doesn't happen? You can download older versions of integron_finder from
  [the releases page](https://github.com/gem-pasteur/integron_finder/releases).
* **Can you reliably reproduce the issue?** If not, provide details about how often the problem happens and under which conditions it normally happens.
* If the problem is related to working with files (e.g. opening and editing files),

### Suggesting Enhancements

This section guides you through submitting an enhancement suggestion for integron_finder,
including completely new features and minor improvements to existing functionality.
Following these guidelines helps maintainers and the community understand your suggestion :pencil:
and find related suggestions :mag_right:.

Before creating enhancement suggestions, please check [this list](#before-submitting-an-enhancement-suggestion)
as you might find out that you don't need to create one.
When you are creating an enhancement suggestion, please [include as many details as possible](#how-do-i-submit-a-good-enhancement-suggestion).
Fill in [the template](.github/ISSUE_TEMPLATE.md), including the steps that you imagine you would take if the feature you're requesting existed.

#### Before Submitting An Enhancement Suggestion

* **Check if you're using [the latest version of integron_finder](https://github.com/gem-pasteur/integron_finder/releases)**.
* **Perform a [cursory search](https://github.com/gem-pasteur/integron_finder/issues?q=is%3Aopen+is%3Aissue+label%3Aenhancement)**
  to see if the enhancement has already been suggested.
  If it has, add a comment to the existing issue instead of opening a new one.

#### How Do I Submit A (Good) Enhancement Suggestion?

Enhancement suggestions are tracked as [GitHub issues](https://guides.github.com/features/issues/).

* **Use a clear and descriptive title** for the issue to identify the suggestion.
* **Provide a step-by-step description of the suggested enhancement** in as many details as possible.
* **Provide specific examples to demonstrate the steps**.
  Include copy/pasteable snippets which you use in those examples, as [Markdown code blocks](https://help.github.com/articles/markdown-basics/#multiple-lines).
* **Describe the current behavior** and **explain which behavior you expected to see instead** and why.
* **Explain why this enhancement would be useful** to most integron_finder users.
* **Specify which version of integron_finder you're using.** You can get the exact version by running `integron_finder --version` in your terminal.
* **Specify the name and version of the OS you're using.**

### Your First Code Contribution

Unsure where to begin contributing to integron_finder? You can start by looking through these `beginner` and `help-wanted` issues:

* [Beginner issues][beginner] - issues which should only require a few lines of code, and a test or two.
* [Help wanted issues][help-wanted] - issues which should be a bit more involved than `beginner` issues.

### Pull Requests

#### Process

1. Follow the existing code style precedent. This does not need to be strictly
   defined as there are many thousands of lines of examples. Note the lack
   of tabs anywhere in the project, parentheses and spacing, curly bracket
   locations, source code layout, variable scoping, etc. and follow the
   project's standards.
2. For any new functionality, please write a test to be added to Continuous
   Integration (Travis) to test it (tests can be found in the `tests/`
   directory).
3. The project's default copyright and header have been included in any new
   source files.
4. Make sure you have implemented a tests and all tests `python tests/run_tests.py`
   succeed before submitting the PR.
5. Is the code human understandable? This can be accomplished via a clear code
   style as well as documentation and/or comments.
6. The pull request will be reviewed by others, and the final merge must be
   done by the Integron_finder project lead.
7. Documentation must be provided if necessary ([next section](#documentation-style-guide))
8. Fill in [the required template](PULL_REQUEST_TEMPLATE.md)
9. Do not include issue numbers in the PR title

### Style guides

### Git Commit Messages

* Limit the first line to 72 characters or less
* Reference issues and pull requests liberally after the first line
* Consider starting the commit message with an applicable emoji:
    * :art: `:art:` when improving the format/structure of the code
    * :sparkles: `:sparkles:` when implement a new feature
    * :memo: `:pencil2:` when writing docs
    * :bug: `:bug:` when fixing a bug
    * :white_check_mark: `:white_check_mark:` when adding tests
    * :fire: `:fire:` when removing code or files
    * :arrow_up: `:arrow_up:` when upgrading dependencies
    * :arrow_down: `:arrow_down:` when downgrading dependencies


### Python Style guide

*integron_finder* (from version 2.0) is written in Python3.
We try to follow the [zen of python](https://www.python.org/dev/peps/pep-0020/) principles,
and we adopted the [pep8](https://www.python.org/dev/peps/pep-0008/) coding style.
Add docstring compatible with sphinx in [restructuredtext]() format.
In docstring describe all parameters and types as well as the return, and error raise
 directly by the method when applicable.

```python
def foo(param_1, param_2):
    """
    foo do bla bla ...

    :param param_1: what is the param 1
    :type param_1: str
    :param str param_2: what is the param 2,
                      when it is a basic type, it can be specified "on line"
    :returns: what is returned by foo
    :rtype: list of strings
    :raise RuntimeError: when param_1 is longer than param_2
    """
    if len(param_1) > len(param_2):
        raise RuntimeError()
    ...
```
INtegron_Finder project adopt ruff as linter, all the code must pass the `ruff check` step.
To avoid to push incorrect code *IF* *use pre-commit*. So perform `pre-commit install` once you set up your environment
(see developper_guide/installation for more details)

### Documentation Style guide

* We use [sphinx](http://www.sphinx-doc.org/en/stable/) with reStructuredText syntax to document the project.
* We have separate user and developer documentation.
* In developer documentation before diving head first in detailed documentation,
  you should try to get an overview of how your code works.
* We have a new page for each new scripts.
* In developer documentation don't forget the api reference section,
  document every classes, methods (even private) or functions.
* In user documentation add examples, command lines, input and output.
