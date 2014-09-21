#!/usr/bin/env python

""" I/O utilities for dealing with SIU files """

import os
import struct
import numpy as np

__all__ = ["read_line_list", "write_line_list"]


def chunk(l, n):
    return [l[i:i + n] for i in range(0, len(l), n)]

def guess_mode(filename):
    textchars = ''.join(map(chr, [7,8,9,10,12,13,27] + range(0x20, 0x100)))
    is_binary_string = lambda bytes: bool(bytes.translate(None, textchars))
    with open(filename) as fp:
        is_binary = is_binary_string(fp.read(1024))
    return ["w", "b"][is_binary]

def _tuple_to_recarray(transition):
    columns = ("wavelength", "excitation_potential", "J_low", "J_up", "xx1",
        "log(gf)", "xx2", "log(c4)?", "log(c6)", "xx3", "xx4", "xx5", "xx6",
        "xx7", "xx8", "xx9", "xx10", "atomic_number", "xx11", "ionised")
    return np.core.records.fromrecords([transition], names=columns, formats=None)
    

def read_line_list(filename, mode=None, ignore_blanks=True):
    """
    Read a line list and return a record array.

    :param filename:
        The path of the line list.

    :type filename:
        str

    :param mode: [optional]
        The mode to open the file in. Specifying mode "w" indicates ASCII format.
        Specifying the mode as "b" indicates binary format. If no mode is given,
        it will be guessed from the file contents.

    :type mode:
        str

    :returns:
        A record array of the line list.
    """

    # Read it into tuples, regardless of binary/ascii
    if mode is None:
        mode = guess_mode(filename)

    if "b" in mode.lower():
        records = _read_binary_line_list(filename)
    else:
        records = _read_ascii_line_list(filename)

    # Convert to a record array
    columns = ("wavelength", "excitation_potential", "J_low", "J_up", "xx1",
        "log(gf)", "xx2", "log(c4)?", "log(c6)", "xx3", "xx4", "xx5", "xx6",
        "xx7", "xx8", "xx9", "xx10", "atomic_number", "xx11", "ionised")
    data = np.core.records.fromrecords(records, names=columns, formats=None)
    return data[data["wavelength"] > 0] if ignore_blanks else data


def write_line_list(filename, data, mode="w", clobber=False):
    """
    Write the line list to disk.

    :param filename:
        The path of the line list to write to.

    :type filename:
        str

    :param data:
        The line list information as a record array.

    :type data:
        :class:`np.core.records.recarray`

    :param mode: [optional] 
        The mode to open the file in. The default (w) indicates ASCII format.
        Specifying the mode as "b" indicates binary format.

    :type mode:
        str

    :param clobber: [optional]
        Clobber the file if it already exists.

    :type clobber:
        bool
    """

    if os.path.exists(filename) and not clobber:
        raise IOError("filename {0} exists and not clobbering".format(filename))

    # Data should be in record array. Can pass to either binary or ascii writer.
    if "b" in mode.lower():
        return _write_binary_line_list(filename, data)
    else:
        return _write_ascii_line_list(filename, data)


def _write_ascii_line_list(filename, data):
    """
    Write line list to ASCII format.
    """

    #  2010.974   0.0300  2.0  3.0   0.1  -3.3800   0.0000   0.0000 -31.1250   0
    #.0000   0.0000         2    uv6   3P   1G  KP   STD   14  0  1
    record = "{0:10.3f}{1:9.4f}{2:5.1f}{3:5.1f}{4:6.1f}{5:9.4f}{6:9.4f}{7:9.4f}"\
        "{8:9.4f}{9:9.4f}{10:9.4f}{11:10.0f}  {12:>5s}  {13:>3s}  {14:>3s} {15"\
        ":>3s}   {16:>3s}  {17:3.0f}{18:3.0f}{19:3.0f}\n"

    scales = np.ones(len(data.dtype.names))
    scales[6] = 1e+8
    
    with open(filename, "w+") as fp:
        for row in data:
            row_data = row * scales
            fp.write(record.format(*row_data))

    return True


def _write_binary_line_list(filename, data, num=200):
    """
    Write line list to binary format.
    """

    assert 3.8e5 >= len(data), "SI line list format supports at most 380000 lines"
    
    num_structures = int(np.ceil(len(data)/float(num)))

    # Create the header
    header_values = np.zeros((1900))
    header_values[:num_structures] = data["wavelength"][::num]
    header = struct.pack("<1900d", *header_values)

    # Create the structures
    contents = header
    for i in xrange(num_structures):
        contents += _pack_binary_structure(data[i*num:(i+1)*num], num)

    with open(filename, "wb+") as fp:
        fp.write(contents)
    return True


def _pack_binary_structure(data, num=200):

    assert num >= len(data)

    columns = ["wavelength", "excitation_potential", "J_low", "J_up", "xx1",
        "log(gf)", "xx2", "log(c4)?", "log(c6)", "xx3", "xx4", "xx5", "xx6",
        "xx7", "xx8", "xx9", "xx10", "atomic_number", "xx11", "ionised"]

    data = data[columns]

    if num > len(data):
        # We will need to copy the data and pad it.
        data = data.copy()
        padded_row = np.array([(0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0, 
            ' ', ' ', ' ', ' ', ' ', 0, 0, 0)], dtype=data.dtype)
        for i in xrange(num - len(data)):
            data = np.append(data, padded_row)

    colnames = data.dtype.names[:12]
    scales = np.ones(len(colnames))
    scales[6] = 1e+8
    formatted_data = struct.pack(
        "<{0}d{0}d{0}f{0}f{0}f{0}f{0}f{0}f{0}f{0}f{0}f{0}l".format(num),
        *np.array([data[c]*s for c, s in zip(colnames, scales)]).flatten())

    sizes = (5, 3, 3, 3, 3)
    colnames = data.dtype.names[12:17]
    formatted_byte_data = struct.pack("<{0}s{1}s{2}s{3}s{4}s".format(
        *[num*s for s in sizes]), *["".join(list(data[c])) for c in colnames])

    colnames = data.dtype.names[17:]
    formatted_element_data = struct.pack("<{0}b{0}b{0}b".format(num),
        *np.array([data[c] for c in colnames]).flatten())

    return "".join([formatted_data, formatted_byte_data, formatted_element_data])


def _read_ascii_line_list(filename):
    """
    Read an ASCII-formatted line list and return a list containing tuples for 
    each record.
    """

    with open(filename, "r") as fp:
        lines = fp.readlines()
    return map(_read_ascii_line, lines)


def _read_ascii_line(line):
    """
    Read the structure containined in an ASCII line and return the record as a
    tuple.
    """

    #  2010.974   0.0300  2.0  3.0   0.1  -3.3800   0.0000   0.0000 -31.1250   \
    #0.0000   0.0000         2    uv6   3P   1G  KP   STD   14  0  1
    record = []
    record.extend(map(float, line.split()[:12]))
    record[6] = record[6] * 1e-8
    record[-1] = int(record[-1])

    # Split the text:
    record.extend(map(str.strip, chunk(line[101:126], 5)))

    # And the final element information
    record.extend(map(int, line[127:].split()))
    return tuple(record)


def _read_binary_line_list(filename):
    """
    Read a binary-formatted line list and return a list containing tuples for
    each record.
    """

    with open(filename, "rb") as fp:
        # Skip the header (this is also length 15200)
        fp.seek(0x3b60)
        contents = fp.read()

    size = 15200
    return sum([_read_binary_structure(contents[i*size:]) \
        for i in xrange(len(contents)/size)], [])
    

def _read_binary_structure(contents, num=200):
    """
    Read the structure contained in a binary line list and return the records
    of the structure.
    """

    data_fmt = "<{0}d{0}d{0}f{0}f{0}f{0}f{0}f{0}f{0}f{0}f{0}f{0}l".format(num)
    data = struct.unpack_from(data_fmt, contents)

    sizes = (5, 3, 3, 3, 3)
    offset = struct.calcsize(data_fmt)
    byte_fmt = "<{0}s{1}s{2}s{3}s{4}s".format(*[num*s for s in sizes])
    str_data = [chunk(d, s) for d, s in \
        zip(struct.unpack_from(byte_fmt, contents[offset:]), sizes)]

    offset += struct.calcsize(byte_fmt)
    element_fmt = "<{0}b{0}b{0}b".format(num)
    element_data = struct.unpack_from(element_fmt, contents[offset:])

    # Join records together
    records = []
    for i in xrange(num):
        record = list(data[i::num])
        record[6] *= 1e-8
        record.extend([each[i] for each in str_data])
        record.extend(element_data[i::num])
        records.append(tuple(record))
    return records
