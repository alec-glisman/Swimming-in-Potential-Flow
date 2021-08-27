# @AUTHOR: alec-glisman (GitHub)
#
# NOTE:
#   This module is built to serve as a style guide for generating
#   publication quality plots.
#   Guiding principles for plot design:
#      * Axes as minimal as possible
#      * Key elements "pop out" more with thicker lines
#      * Minimize negative space
#
# USEAGE:
#    # Instantiate class instance anc create back-end figure
#    sample_plot = PlotStyling(num_lines,
#                        r"x_axis_label", r"y_axis_label",
#                        title=None, loglog=False,
#                        outputDir="", figName="sample_plot",
#                        eps=False,
#                        continuousColors=False)
#    mod_plot.make_plot()
#
#    # Curves & scatters
#    mod_plot.curve(x_data, y_data, thin_curve=False, dashed_curve=False,
#                zorder=0, label=r"$data$", color=None)
#
#    # Legend & text elements
#    mod_plot.legend(title=None, loc='best',
#                    bbox_to_anchor=(0.01, 0.01, 0.98, 0.98),
#                    ncol=1)
#
#    # Save plot
#    mod_plot.save_plot()


# SECTION: Dependencies

# System level information
import sys                             # reading and redirecting std channels
import os                              # opening dev/null channel

# Warnings
import logging                         # overwriting matplotlib errors
import warnings                        # catch mathTex warnings
from contextlib import contextmanager  # decorator to suppress std out

# Data structures
import numpy as np                     # general tensor data-structures

# Plotting
import matplotlib as mpl               # plt parent library to set global parameters
import matplotlib.pyplot as plt        # plotting library
import matplotlib.ticker as mtk        # matplotlib ticking attributes

# Colorschemes
import colorcet as cc                  # plotting colormaps

# !SECTION


class PlotStyling:

    # SECTION: Static parameters

    # Max size in matplotlib to prevent overflow error
    # @Source: https://stackoverflow.com/questions/37470734/matplotlib-giving-error-overflowerror-in-draw-path-exceeded-cell-block-limit
    mpl.rcParams['agg.path.chunksize'] = 10000000

    # Overall plt parameters (@AUTHOR is colorblind)
    # either 'tableau-colorblind10' or 'seaborn-colorblind'
    plt.style.use(["tableau-colorblind10"])
    # Set the packages required for LaTeX rendering
    mpl.rcParams['text.latex.preamble'] = r"\usepackage{amsfonts,amsmath,amsbsy,amssymb,bm,amsthm,mathrsfs,fixmath}"

    # Use mathtext, not LaTeX
    mpl.rcParams.update({
        'text.usetex': True,

        # Use the Computer modern font
        'font.family': 'serif',
        'font.serif': 'cmr10',
        'mathtext.fontset': 'cm',

        # Use ASCII minus
        'axes.unicode_minus': False,
    })

    # Set text element font sizes
    title_size = 40
    axis_label_size = 40
    tick_size = 30
    legend_size = 22
    legend_title_size = 26

    # Set padding to use on x and y axis labels
    label_pad = 10

    # Set plotting element sizes
    line_width = 3
    thin_line_width = 2
    marker_size = 28
    dashes = [5, 2, 1, 2]  # 5 points on, 2 off, 1 on, 2 off

    # Linear color options
    col_diverge = cc.bky
    col_sunset = cc.CET_L8
    col_beetle = mpl.cm.get_cmap('viridis').colors

    # Categorical colors
    darkBlue = cc.glasbey_dark[10]
    black = 'k'
    red = cc.glasbey_dark[0]
    darkgreen = cc.glasbey_dark[2]
    magenta = cc.glasbey_dark[13]
    orange = cc.glasbey_dark[41]
    lightBlue = cc.glasbey_dark[16]
    lightGreen = cc.glasbey_dark[44]
    purple = cc.glasbey_dark[42]
    brown = cc.glasbey_dark[15]
    col_categorical = [darkBlue, red, darkgreen, magenta, black,
                       orange, lightBlue, lightGreen, purple, brown]


# !SECTION


# SECTION: Initializer

    def __init__(self, numLines,
                 x_label, y_label,
                 title=None, loglog=False,
                 outputDir=None, figName=None, eps=False,
                 continuousColors=False):

        # Plot text elements
        self.x_label = x_label
        self.y_label = y_label
        self.title = title

        # Overall figure formatting
        self.numLines = None  # how many data-series will be plotted on current figure
        self.loglog = loglog  # Whether to output log-log scaled axes

        # Handles for the plot (from matplotlib)
        self.fig = None
        self.ax = None

        # Plot I/O
        self.figName = figName      # filename of plot
        self.outputDir = outputDir  # directory to save plot
        self.eps = eps              # whether to output in .eps or .png format

        # Set color scheme
        # REVIEW: There are 3 options for "continuous" colorshemes
        #           * Diverging colors: 'cc.bky'                             (blue-black-brown)
        #           * Linear colors 1:  'cc.CET_L8'                          (blue-purple-pink-orange-yellow)
        #           * Linear colors 2:  'mpl.cm.get_cmap('viridis').colors'  (purple-teal-green-yellow)

        # Colorscheme
        self.continuous_colors = continuousColors
        self.color_ID = 0
        self.color_ID_increment = 1
        self.color_map = None

        if self.continuous_colors:

            # Input checking
            if (numLines > 255):
                raise ValueError(
                    "numLines input must be < 255 for continuous colorscheme output")

            # Choose a map
            self.color_map = self.col_sunset
            self.color_ID_increment = int(255 / numLines)

        else:

            # Input checking
            if (numLines > len(self.col_categorical)):
                raise ValueError(
                    f"numLines input must be < {len(self.col_categorical)} for continuous colorscheme output")

            # Choose a map
            self.color_map = self.col_categorical

# !SECTION


# SECTION: Make the plot figure

    def make_plot(self, showPlot=False):

        if showPlot:

            plt.switch_backend('module://ipykernel.pylab.backend_inline')
            # Lower-resolution figure
            figsize = (10, 7)
            dpi = 400

        else:

            # Does not actually open any GUI so that it can be run on servers
            plt.switch_backend('Agg')
            # High-resolution figure
            figsize = (10, 7)
            dpi = 600

        # Construct plot
        self.fig, self.ax = plt.subplots(figsize=figsize, dpi=dpi)

        # Set axis labels
        self.ax.set_xlabel(
            self.x_label, fontsize=self.axis_label_size, labelpad=self.label_pad)
        self.ax.set_ylabel(
            self.y_label, fontsize=self.axis_label_size, labelpad=self.label_pad)
        # Set title
        if self.title is not None:
            self.ax.set_title(
                self.title, fontsize=self.title_size, fontweight='bold')

        # Update axis ticks
        self.ax.xaxis.set_tick_params(labelsize=self.tick_size, which="major",
                                      direction='in', width=2, length=7)
        self.ax.xaxis.set_tick_params(labelsize=self.tick_size, which="minor",
                                      direction='in', width=1.5, length=7)
        self.ax.yaxis.set_tick_params(labelsize=self.tick_size, which="major",
                                      direction='in', width=2, length=7)
        self.ax.yaxis.set_tick_params(labelsize=self.tick_size, which="minor",
                                      direction='in', width=1.5, length=7)
        self.ax.ticklabel_format(useMathText=True)
        self.ax.yaxis.set_ticks_position('both')
        self.ax.xaxis.set_ticks_position('both')

        # Add minor ticks at all integers for log-log plots
        if self.loglog:
            pass

        # Make border larger
        self.ax.spines["top"].set_linewidth(2)
        self.ax.spines["bottom"].set_linewidth(2)
        self.ax.spines["left"].set_linewidth(2)
        self.ax.spines["right"].set_linewidth(2)

# !SECTION


# SECTION: Continuous curves

    def curve(self, x_data, y_data,
              thin_curve=False, dashed_curve=False,
              zorder=None, label=None, color=None):

        # Set color
        local_color = color
        if local_color is None:
            local_color = self.color_map[self.color_ID]
            self.color_ID += self.color_ID_increment        # Increment color ID

        # Set line style
        local_curve = '-'  # solid line style
        if dashed_curve:
            local_curve = '--'  # dashed line style

        # Set linewidth
        local_linewidth = self.line_width
        if thin_curve:
            local_linewidth = self.thin_line_width

        # Plot input data
        if self.loglog:
            self.ax.loglog(x_data, y_data, local_curve,
                           color=local_color, label=label,
                           linewidth=local_linewidth, zorder=zorder)
        else:
            self.ax.plot(x_data, y_data, local_curve,
                         color=local_color, label=label,
                         linewidth=local_linewidth, zorder=zorder)

        self.fix_unicode_chars()  # Fix any character rendering errors

# !SECTION


# SECTION: Discrete data points

    def scatter(self, x_data, y_data,
                marker=None, zorder=None, label=None, color=None):

        # Set color
        local_color = color
        if local_color is None:
            local_color = self.color_map[self.color_ID]
            self.color_ID += self.color_ID_increment        # Increment color ID

        # Set marker
        local_marker = marker
        if local_marker is None:
            local_marker = 'o'

        # Plot input data
        if self.loglog:
            scatter = self.ax.loglog(x_data, y_data, local_marker,
                                     color=local_color, label=label,
                                     linewidth=self.line_width, zorder=zorder)
        else:
            scatter = self.ax.plot(x_data, y_data, local_marker,
                                   color=local_color, label=label,
                                   linewidth=self.line_width, zorder=zorder)

        self.fix_unicode_chars()  # Fix any character rendering errors

        return scatter

    def scatter_dashed(self, x_data, y_data, marker=None, zorder=None, label=None, color=None):

        # Set color
        local_color = color
        if local_color is None:
            local_color = self.color_map[self.color_ID]
            self.color_ID += self.color_ID_increment        # Increment color ID

        # Set marker
        local_marker = marker
        if local_marker is None:
            local_marker = 'o'

        # Plot input data
        if self.loglog:
            scatter, = self.ax.loglog(x_data, y_data, local_marker,
                                      color=local_color, label=label,
                                      linewidth=self.line_width, zorder=zorder)
        else:
            scatter, = self.ax.plot(x_data, y_data, local_marker,
                                    color=local_color, label=label,
                                    linewidth=self.line_width, zorder=zorder)

        scatter.set_dashes(self.dashes)  # Add connecting dashes

        self.fix_unicode_chars()  # Fix any character rendering errors

        return scatter

# !SECTION


# SECTION: Legend

    def legend(self, title=None, loc='best', bbox_to_anchor=None, ncol=1):

        legend = self.ax.legend(title=title, loc=loc, bbox_to_anchor=bbox_to_anchor,
                                fontsize=self.legend_size,
                                frameon=False, fancybox=True,
                                framealpha=0.90, shadow=False,
                                edgecolor='k',
                                ncol=ncol)

        legend.get_title().set_fontsize(self.legend_title_size)

# !SECTION


# SECTION: Text elements


    def textbox(self, text, x_loc, y_loc,
                horz_align="center", vert_align="center"):

        plt.text(x_loc, y_loc, text,
                 fontsize=self.legend_title_size,
                 ha=horz_align, va=vert_align)

# SECTION


# SECTION: Axes

    def set_yaxis_scientific(self):

        plt.ticklabel_format(style='sci', axis='y',
                             scilimits=(0, 2), useOffset=True)
        t = self.ax.yaxis.get_offset_text()
        t.set_size(self.tick_size)

    def set_xaxis_scientific(self):

        plt.ticklabel_format(style='sci', axis='x',
                             scilimits=(0, 2), useOffset=True)
        t = self.ax.xaxis.get_offset_text()
        t.set_size(self.tick_size)

    def set_xaxis_tick_scalar_formatter(self):

        self.ax.xaxis.set_minor_formatter(mtk.ScalarFormatter())
        self.ax.xaxis.set_major_formatter(mtk.ScalarFormatter())

    def set_yaxis_tick_scalar_formatter(self):

        self.ax.yaxis.set_minor_formatter(mtk.ScalarFormatter())
        self.ax.yaxis.set_major_formatter(mtk.ScalarFormatter())

    def set_major_minor_ticks(self, xMajorLoc=None, xMinorLoc=None, yMajorLoc=None, yMinorLoc=None):

        if xMajorLoc is not None:
            self.ax.xaxis.set_major_locator(
                mpl.ticker.MultipleLocator(xMajorLoc))

        if xMinorLoc is not None:
            self.ax.xaxis.set_minor_locator(
                mpl.ticker.MultipleLocator(xMinorLoc))

        if yMajorLoc is not None:
            self.ax.yaxis.set_major_locator(
                mpl.ticker.MultipleLocator(yMajorLoc))

        if yMinorLoc is not None:
            self.ax.yaxis.set_minor_locator(
                mpl.ticker.MultipleLocator(yMinorLoc))

    def set_latex_axis_labels(self, axis='x', axisScale=np.pi, denominator=2, minorTicks=4, latex=r'\pi'):

        major = Multiple(denominator, axisScale, latex)
        minor = Multiple(denominator*minorTicks, axisScale, latex)

        # Output results to figure
        if (axis == 'x'):
            self.ax.xaxis.set_major_locator(major.locator())
            self.ax.xaxis.set_minor_locator(minor.locator())
            self.ax.xaxis.set_major_formatter(major.formatter())

        elif (axis == 'y'):
            self.ax.yaxis.set_major_locator(major.locator())
            self.ax.yaxis.set_minor_locator(minor.locator())
            self.ax.yaxis.set_major_formatter(major.formatter())

# !SECTION


# SECTION: Save figure

    def save_plot(self, showPlot=False):

        self.fix_unicode_chars()

        if (self.figName is not None) and (self.eps == True):
            with self.suppress_stdout():
                plt.savefig(self.outputDir + self.figName + '.eps',
                            format='eps', bbox_inches='tight')

        elif (self.figName is not None) and (self.eps == False):
            plt.savefig(self.outputDir + self.figName + '.png',
                        format='png', bbox_inches='tight')

        if showPlot:
            plt.show()

# !SECTION


# SECTION: Remove unicode issues

    # Removes issue with \times character in mathTex output of plots:
    # @Source: https://stackoverflow.com/questions/47253462/matplotlib-2-mathtext-glyph-errors-in-tick-labels

    def fix_unicode_chars(self):

        # Force the figure to be drawn
        logger = logging.getLogger('matplotlib.mathtext')
        original_level = logger.getEffectiveLevel()
        logger.setLevel(logging.ERROR)

        with warnings.catch_warnings():

            warnings.simplefilter(
                'ignore', category=mpl.mathtext.MathTextWarning)
            self.fig.canvas.draw()

        logger.setLevel(original_level)

        # Remove '\mathdefault' from all minor tick labels
        # 'with' block added by @Alec in December 2020 to avoid UserWarning related to relabling axis
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore")

            labels = [label.get_text().replace(r'\mathdefault', '')
                      for label in self.ax.get_xminorticklabels()]
            self.ax.set_xticklabels(labels, minor=True)

            labels = [label.get_text().replace(r'\mathdefault', '')
                      for label in self.ax.get_yminorticklabels()]
            self.ax.set_yticklabels(labels, minor=True)

        # Added by @Alec
        self.ax.xaxis.set_tick_params(
            which='both', labelsize=self.tick_size, pad=self.label_pad)
        self.ax.yaxis.set_tick_params(
            which='both', labelsize=self.tick_size, pad=self.label_pad)

    # Suppresses console output of function calls (used for eps no transparent background issue)
    # @Source: https://stackoverflow.com/questions/2125702/how-to-suppress-console-output-in-python
    # @User: charleslparker
    @ contextmanager
    def suppress_stdout(self):

        with open(os.devnull, "w") as devnull:
            old_stdout = sys.stdout
            sys.stdout = devnull
            old_stderr = sys.stderr
            sys.stderr = devnull

            try:
                yield

            finally:
                sys.stdout = old_stdout
                sys.stderr = old_stderr

# !SECTION


# SECTION: Get \pi in axis labels of plots

# @Source: https://stackoverflow.com/questions/40642061/how-to-set-axis-ticks-in-multiples-of-pi-python-matplotlib
# @User: Scott Centoni
def multiple_formatter(denominator=2, number=np.pi, latex=r'\pi'):

    def gcd(a, b):

        while b:
            a, b = b, a % b

        return a

    def _multiple_formatter(x, pos):

        den = denominator
        num = np.int(np.rint(den*x/number))
        com = gcd(num, den)
        (num, den) = (int(num/com), int(den/com))

        if den == 1:

            if num == 0:
                return r'$0$'

            if num == 1:
                return r'$%s$' % latex

            elif num == -1:
                return r'$-%s$' % latex

            else:
                return r'$%s%s$' % (num, latex)
        else:

            if num == 1:
                return r'$\frac{%s}{%s}$' % (latex, den)

            elif num == -1:
                return r'$\frac{-%s}{%s}$' % (latex, den)

            else:
                return r'$\frac{%s%s}{%s}$' % (num, latex, den)

    return _multiple_formatter


class Multiple:

    def __init__(self, denominator=2, number=np.pi, latex=r'\pi'):

        self.denominator = denominator
        self.number = number
        self.latex = latex

    def locator(self):
        return plt.MultipleLocator(self.number / self.denominator)

    def formatter(self):
        return plt.FuncFormatter(multiple_formatter(self.denominator, self.number, self.latex))

# !SECTION
