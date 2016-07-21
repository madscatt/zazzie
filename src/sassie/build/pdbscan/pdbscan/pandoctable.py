# -*- coding: utf-8 -*-
"""
Pandoctable:

Methods for creating tables in Pandoc Markdown
"""

'''
    SASSIE: Copyright (C) 2011 Joseph E. Curtis, Ph.D.

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
'''

from textwrap import TextWrapper


def create_pandoc_table(header, contents, widths, just):
    """
    Split text for several columns to provide list of lists containing strings
    of the correct width for the colmn in a table.

    @type  header:    list
    @param header:    List of strings (one header per column)
    @type  contents:  list
    @param contents:  List of lists (sub-list contains the data for each column
                      in each row)
    @type  widths:    list
    @param widths:    Widths for each column
    @type  just:      list
    @param just:      List of justification choices [(l)eft, (r)ight, (c)entre]
                      for each column (one element per column)
    @rtype:           string
    @return:          Entire table with lines separated by '\\n' character
    """

    table = []

    if header:
        full_width = sum(widths) + len(widths) - 1
        table.append('-' * full_width)
        for line in create_table_lines(header, widths, just, isheader=True):
            table.append(line)

    col_indicators = create_table_col_line(widths)

    table.append(col_indicators)

    last = len(contents)

    for i in range(0, last):

        lines = create_table_lines(contents[i], widths, just)

        if i == last - 1:

            if not lines[-1].strip():

                # Blank final lines must not be printed
                lines = lines[:-1]

        for line in lines:
            # if not line.strip(): continue ## @NOTE: ZHL added to skip empty
            # lines for compact report
            table.append(line)

    table.append(col_indicators)
    table.append('\n')

    return table


def wrap_columns(cols, widths):
    """
    Split text for several columns to provide list of lists containing strings
    of the correct width for the colmn in a table.

    @type  cols:      list
    @param cols:      List of strings (one element per column)
    @type  widths:    list
    @param widths:    Widths for each column
    @rtype:           list
    @return:          list of lists containing lines wrapped to the width of
                      each column

    """

    wrap_cols = []

    for i in range(0, len(widths)):

        if len(cols[i]) > widths[i]:

            wrapper = TextWrapper(width=widths[i])
            wrap_cols.append(wrapper.wrap(cols[i]))

        else:

            wrap_cols.append([cols[i]])

    return wrap_cols


def create_table_lines(cols, widths, just, isheader=False):
    """
    Prepare text for a set of columns for output as a Pandoc Markdown multiline
    table. This means that rows can span multiple lines but must be separated
    by a blank line. For header lines no trailing space is added. Output is as
    list of lines where columns are separated by a single space. Strings are
    padded with spaces if they do not fill the width allocated to them.

    @type  cols:        list
    @param cols:        List of strings (one element per column)
    @type  widths:      list
    @param widths:      Widths for each column
    @type  just:        list
    @param just:        Justification for each column 'l'/'r'/'c' for left,
                        right or centre
    @type  isheader:    boolean
    @param isheader:    Indicates if the goal is to generate a header line
    @rtype:             list
    @return:            list of strings representing the lines making up a row
                        of a Pandoc Markdown multiline table
    """

    # Wrap each column to fit assigned width
    wrap_cols = wrap_columns(cols, widths)

    no_rows = max(len(x) for x in wrap_cols)

    for i in range(0, len(wrap_cols)):

        width = widths[i]

        # Justify each line
        for j in range(0, len(wrap_cols[i])):

            if just[i] == 'r':

                wrap_cols[i][j] = wrap_cols[i][j].rjust(width)

            elif just[i] == 'c':

                wrap_cols[i][j] = wrap_cols[i][j].center(width)

            else:

                wrap_cols[i][j] = wrap_cols[i][j].ljust(width)

        # Pad line to fit column if necessary
        wrap_cols[i] += [" " * width] * (no_rows - len(wrap_cols[i]))

    table_lines = []

    for row in zip(*wrap_cols):
        table_lines.append(' '.join(row))

    if not isheader:
        table_lines.append('')

    return table_lines


def create_table_col_line(widths):
    """
    Create a dashed header line to show the locations of the columns according
    to the Pandoc Markdown table syntax. Example with columns of
    widths = [5,7,3]:
    ----- ------- ---

    @type  widths:      list
    @param widths:      Widths for each column
    @rtype:             string
    @return:            Line of dashes designating column widths in a Pandoc
                        Markdown standard table
    """

    header = ''

    for width in widths:
        col_head = '-' * width + ' '
        header += col_head

    return header
