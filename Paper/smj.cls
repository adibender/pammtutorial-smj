%%
%% This is file `smj.cls',
\def\fileversion{0.2}
\def\filename{smj}
\def\filedate{2014/02/14}
%%
%% Package `smj' to use with LaTeX2e for preparation of manuscripts
%% to be submitted to Statistical Modelling (http://smj.sagepub.com)
%% License: GPL-3
%% Copyright: (C) Arnost Komarek
%% Please report errors to smj-komarek@karlin.mff.cuni.cz
%%
%% 0.2:  support for bibTeX added
%%
\NeedsTeXFormat{LaTeX2e}
\ProvidesClass{smj}[\filedate\space\fileversion\space\filename\space class by Arnost Komarek]


%%%%% Options
%%%%% ------------------------------------------------------------
\newif\if@submit
\newif\if@publish
\newif\if@shortnames

\@submittrue
\@publishfalse
\@shortnamestrue

\DeclareOption{submit}{
  \@submittrue
  \@publishfalse
}

\DeclareOption{publish}{
  \@submitfalse
  \@publishtrue
}

\DeclareOption{shortnames}{
  \@shortnamestrue
}


%%%%% Process options
%%%%% -----------------------------------------------------------
\ProcessOptions

%% Basic class and its options
\if@submit
  \LoadClass[12pt, a4paper, twoside]{article}
\fi
\if@publish
  \LoadClass[11pt, a4paper, twoside]{article}
\fi

%% Required packages
\RequirePackage{graphicx}
\RequirePackage{color}
\RequirePackage{fancyvrb}
%\RequirePackage{footmisc}
\RequirePackage{amsmath,amssymb}
\RequirePackage{calc}
%%\RequirePackage[T1]{fontenc}
\RequirePackage{url}
\IfFileExists{upquote.sty}{\RequirePackage{upquote}}{}
\definecolor{red2}{rgb}{0.9333333, 0, 0}
\definecolor{mediumblue}{rgb}{0, 0, 0.8039216}
\definecolor{magenta4}{rgb}{0.545098, 0, 0.545098}
\RequirePackage[colorlinks, linkcolor=red2, citecolor=mediumblue, urlcolor=magenta4, breaklinks=true]{hyperref}

\if@submit
  \setlength{\textwidth}{6in}
  \setlength{\textheight}{9in}
  \setlength{\oddsidemargin}{(\paperwidth - \textwidth)/2 - 1in}
  \setlength{\evensidemargin}{(\paperwidth - \textwidth)/2 - 1in}
  \setlength{\voffset}{(\paperheight - \textheight)/2 - 2in}
  \setlength{\parindent}{0pt}
  \setlength{\parskip}{1em}

  \renewcommand{\baselinestretch}{1.8}

  %\RequirePackage{lineno}
  %\renewcommand\thelinenumber{\scriptsize\color{blue}\arabic{linenumber}\hspace{2em}}
  %\linenumbers
\fi
\if@publish
  \setlength{\textwidth}{6in}
  \setlength{\textheight}{8.5in}
  \setlength{\oddsidemargin}{(\paperwidth - \textwidth)/2 - 1in}
  \setlength{\evensidemargin}{(\paperwidth - \textwidth)/2 - 1in}
  \setlength{\voffset}{(\paperheight - \textheight)/2 - 1.5in}
  \setlength{\parindent}{0pt}
  \setlength{\parskip}{1em}
\fi


%% Bibliography
\if@shortnames
  \usepackage[authoryear,round]{natbib}
\else
  \usepackage[authoryear,round,longnamesfirst]{natbib}
\fi
\bibpunct{(}{)}{;}{a}{,}{,}
\bibliographystyle{smj}


%% Author, title, etc.
\newcommand{\Author}[1]{\def\@Author{#1}}
\newcommand{\AuthorRunning}[1]{\def\@AuthorRunning{#1}}

\newcommand{\Affil}[1]{$^{\;\mathbf{#1}}$}

\newcommand{\Title}[1]{\def\@Title{#1}}
\newcommand{\TitleRunning}[1]{\def\@TitleRunning{#1}}

\newcommand{\CorrAddress}[1]{\def\@CorrAddress{#1}}
\newcommand{\CorrEmail}[1]{\def\@CorrEmail{#1}}
\newcommand{\CorrPhone}[1]{\def\@CorrPhone{#1}}
\newcommand{\CorrFax}[1]{\def\@CorrFax{#1}}

\newcommand{\Abstract}[1]{\def\@Abstract{#1}}
\newcommand{\Keywords}[1]{\def\@Keywords{#1}}

\newcommand{\Year}[1]{\def\@Year{#1}}
\newcommand{\Volume}[1]{\def\@Volume{#1}}
\newcommand{\Issue}[1]{\def\@Issue{#1}}
\newcommand{\DOI}[1]{\def\@DOI{#1}}

\newcommand{\SubmitDate}[1]{\def\@SubmitDate{#1}}
\newcommand{\AcceptDate}[1]{\def\@AcceptDate{#1}}

\newcounter{AffilCount}              % Counter for affiliations
\newenvironment{AffilList}{          % A  list environment for the affiliations
  \begin{list}
  {$^{\arabic{AffilCount}}$}
  {\usecounter{AffilCount}
   \setlength{\topsep}{0mm}
   \setlength{\itemsep}{0mm}
   \setlength{\parskip}{0pt}
   \setlength{\parsep}{0pt}
   \setlength{\labelsep}{5pt}
   \setlength{\labelwidth}{1cm}
   \setlength{\leftmargin}{0.3cm} 
   \rmfamily
  }
}{
  \end{list}
}                 % end environment{AfEnum}

\newcommand{\Affiliations}[1]{
  \def\@Affiliations{
    \begin{AffilList}
    #1
    \end{AffilList}
  }
}



%% Defaults
\Author{Firstname Lastname}
\AuthorRunning{Firstname Lastname}

\Title{Title}
\TitleRunning{Running title}

\CorrAddress{Firstname Lastname, Department, University, Street 1, ZIP City, Country}
\CorrEmail{name@university.edu}
\CorrPhone{(+1)$\;$999$\;$888$\;$777}
\CorrFax{(+1)$\;$999$\;$888$\;$666}

\Abstract{+++ !!! +++ Abstracts must be given. +++ !!! +++}
\Keywords{+++ !!! +++ Keywords must be given. +++ !!! +++}

\Year{2999}
\Volume{00}
\Issue{0}
\DOI{10.1177/YYYYYYYXZZZZZZZZ}

\SubmitDate{2999-12-31}
\AcceptDate{2999-12-31}


%% Title page
\renewcommand{\maketitle}{
  \thispagestyle{empty}
  \markboth{\leftline{\hspace*{2em}\@AuthorRunning}}{\rightline{\@TitleRunning\hspace{2em}}}
  \pagestyle{myheadings}

  \begin{flushleft}\sffamily\bfseries
  \LARGE\@Title

  \vspace{0.2ex}
  \Large\@Author
  %\renewcommand{\@makefnmark}{}
  %\setlength{\footnotemargin}{0pt}
  %\setlength{\skip\footins}{5em}
  %\footnote{\small Address for correspondence: \@CorrAddress. E-mail: \@CorrEmail.
  %Phone: \@CorrPhone. Fax: \@CorrFax.}
  %\setcounter{footnote}{0}
  \end{flushleft}
  
  \vspace{0.2ex}
  \@Affiliations

  \vspace{0ex}
  \rule{\textwidth}{1pt}
  \textbf{Address for correspondence:} \@CorrAddress. 

  \vspace{-1em}
  \textbf{E-mail:} \texttt{\@CorrEmail}.

  \vspace{-1em}
  \textbf{Phone:} \@CorrPhone. 

  \vspace{-1em}
  \textbf{Fax:} \@CorrFax.

  \vspace{0ex}
  \rule{\textwidth}{1pt}
  \textbf{Abstract:}\enspace \@Abstract

  \vspace{0ex}
  \rule{\textwidth}{1pt}
  \textbf{Key words:}\enspace \@Keywords

}

%% Sectioning
\renewcommand{\thesection}{\arabic{section}}
\renewcommand{\section}{%
  \@startsection{section}{1}{\z@}%
                {-3.5ex \@plus -1ex \@minus -.2ex}%
                {2.3ex \@plus.2ex}%
                {\sffamily\Large\bfseries}%
}

\renewcommand{\thesubsection}{\arabic{section}.\arabic{subsection}}
\renewcommand{\subsection}{%
  \@startsection{subsection}{2}{\z@}%
                {-0.25ex\@plus -1ex \@minus -.1ex}%
                {1.5ex \@plus .2ex}%
                {\rmfamily\large\bfseries}%
}

\renewcommand{\thesubsubsection}{\arabic{section}.\arabic{subsection}.\arabic{subsubsection}}
\renewcommand{\subsubsection}{%
  \@startsection{subsubsection}{3}{\z@}%
                {-0.1ex\@plus -1ex \@minus -.1ex}%
                {1.5ex \@plus .2ex}%
                {\itshape\normalsize}%
}

%% Numbering of equations
\renewcommand{\theequation}{\thesection.\arabic{equation}}
\numberwithin{equation}{section}


%% Some commands
%%%\providecommand{\url}[1]{\texttt{\textcolor{magenta4}{#1}}}

\endinput
%%
%% End of file `smj.cls'.
