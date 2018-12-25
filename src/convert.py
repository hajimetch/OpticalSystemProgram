# -*- coding: utf-8 -*-
"""
Created on Fri Mar 16 09:59:09 2018

@author: hajimetch
"""

import configparser
from enum import Enum

#
# Associate setting items with each string expression and number.
#


class NumberedEnum(Enum):
    def __new__(cls, value, number):
        obj = object.__new__(cls)
        obj._value_ = value
        obj.number = number
        return obj


class Projection(Enum):
    SPOT = 'spot'
    PATTERN = 'pattern'


class Reflector(Enum):
    TRUE = 'true'
    FALSE = 'false'


class Telephoto(NumberedEnum):
    TRUE = ('true', 1)
    FALSE = ('false', 0)


class InputDataType(NumberedEnum):
    SURFACE = ('surface data', 0)
    POWER = ('power data', 1)
    POWER_LIGHT = ('power and light data', 2)
    GROUP = ('group data', 3)


class TracingNumber(NumberedEnum):
    SINGLE = ('single', 0)
    MULTI = ('multi', 5)


class IntervalAdjustment(NumberedEnum):
    NONE = ('none', 0)
    FOCAL_LENGTH = ('specify focal length', 1)
    FOCAL_POSITION = ('specify focal position', 2)
    TELEPHOTO_POWER_0 = ('telephoto system power set to 0', 3)
    FOCAL_POSITION_LENGTH = ('specify focal position and length', 4)
    FOCAL_POSITION_SEPARATED = ('specify focal position of separated system',
                                4)


class ObjectDistance(NumberedEnum):
    INFINITY = ('to infinity', 0)
    FIRST_POWER = ('to first power', 1)
    FIRST_SURFACE = ('to first surface', 2)
    PRINCIPAL_POSITION = ('to principal position of projection side', 3)


class IncidentLightRange(NumberedEnum):
    SYSTEM_SOURCE = ('by system and light source', 1)
    APERTURE_IRIS_SOURCE = ('by lens aperture iris and light source', 2)
    APERTURE_IRIS = ('by lens aperture and iris', 3)


class TracingDirection(NumberedEnum):
    SCREEN_TO_APERTURE = ('screen to aperture', 0)
    SCREEN_TO_APERTURE_EXT = ('screen to aperture extended', -1)
    SOURCE_TO_APERTURE = ('light source to aperture', 1)
    SOURCE_TO_APERTURE_EXT = ('light source to aperture extended', 2)


class SplitPattern(NumberedEnum):
    NONE = ('none', 0)
    FIRST_INTERVAL = ('first power and interval', 1)
    LAST_INTERVAL = ('last power and interval', 2)
    RATIO_INTERVAL = ('power ratio and interval', 3)
    FIRST_LAST = ('first and last powers', 4)


class PTLCSpecification(NumberedEnum):
    NONE = ('none', -1)
    FRONT = ('front', 0)
    BACK = ('back', 1)
    FRONT_OVER_BACK = ('front/back', 2)
    BACK_OVER_FRONT = ('back/front', 3)


class ScreenDisplayDimension(NumberedEnum):
    XYPLANE_COORDINATES = ('XY-plane and ray passing point coordinates', 0)
    XYPLANE_3D = ('XY-plane and 3D', 1)
    XYPLANE_XZPLANE = ('XY-plane and XZ-plane', 2)
    SPLIT = ('split view', 3)


#
# Define Exceptions.
#


class ElementBoundsException(Exception):
    pass


#
# Define general methods.
#


def write_dataline(fout, config, args):
    data = []
    for e in args:
        data.append(get_data(fout, config, e))
    fout.write(','.join(data) + '\n')


def get_data(fout, config, e):
    if len(e) == 1:
        return e[0]
    elif len(e) == 2:
        return config.get(e[0], e[1], fallback='')
    elif len(e) == 3:
        return config.get(e[0], e[1], fallback=e[2])
    else:
        raise ElementBoundsException("Unexpected number of elements")


def one_of_valid(config, args):
    for e in args:
        if (config.get(e[0], e[1], fallback='') != ''):
            return [e]
    return []


def all_of_valid(config, args):
    ones = []
    for e in args:
        if (config.get(e[0], e[1], fallback='') != ''):
            ones.append(e)
    return ones


#
# Main part of this script.
#

config = configparser.ConfigParser()
config.read('example.ini')

projection = Projection(
    config.get('Parameters', 'projection', fallback='spot'))
reflector = Reflector(config.get('Parameters', 'reflector', fallback='false'))
telephoto = Telephoto(config.get('Parameters', 'telephoto', fallback='false'))
input_data_type = InputDataType(
    config.get('Parameters', 'input data type', fallback='surface data'))
tracing_number = TracingNumber(
    config.get('Parameters', 'tracing number', fallback='single'))
interval_adjustment = IntervalAdjustment(
    config.get('Parameters', 'interval adjustment', fallback='none'))
object_distance = ObjectDistance(
    config.get(
        'Parameters', 'object distance specification', fallback='to infinity'))
incident_light_range = IncidentLightRange(
    config.get(
        'Parameters',
        'incident light range definition',
        fallback='by system and light source'))
tracing_direction = TracingDirection(
    config.get(
        'Parameters', 'tracing direction', fallback='screen to aperture'))
screen_display_dimension = ScreenDisplayDimension(
    config.get(
        'View Parameters',
        'screen display dimension',
        fallback='XY-plane and ray passing point coordinates'))

split_patterns = []
for e in config.get('System', 'split patterns', fallback='none').split(','):
    split_patterns.append(SplitPattern(e.strip()))

ptlc_specifications = []
for e in config.get(
        'System', 'power to lens curvature specifications',
        fallback='none').split(','):
    ptlc_specifications.append(PTLCSpecification(e.strip()))

with open('datafile.txt', 'wt') as fout:
    fout.write('*SDT:data ')
    system_name = 'FS' if projection == Projection.SPOT else ''
    system_name += 'REF' if reflector == Reflector.TRUE else ''
    system_name += config.get('Parameters', 'system name').strip()
    write_dataline(fout, config, [[system_name]])

    fout.write('*IAF:data ')
    write_dataline(fout, config, [[str(telephoto.number)]])

    fout.write('*MSDT:data ')
    args = []
    args.append(['Telephoto', 'distance between systems'])
    args.append(['Telephoto', 'focal length'])
    args.append(['Telephoto', 'principal position of object side'])
    args.append(['Telephoto', 'principal position of projection side'])
    args.append(['Telephoto', 'aperture for display'])
    args.append(['Telephoto', 'refractive index of projection region'])
    args.append(['Telephoto', 'adjusting interval number'])
    write_dataline(fout, config, args)

    fout.write('*SCL:data ')
    write_dataline(fout, config, [['Parameters', 'scale', '100']])

    fout.write('*NPSL:data ')
    args = []
    args.append(['System', 'number of powers'])
    args.append(['System', 'number of surfaces'])
    args.append(['System', 'number of lenses'])
    if (input_data_type == InputDataType.GROUP):
        args.append(['System', 'number of groups'])
    ones = []
    ones.append(['System', 'adjusting power interval number'])
    ones.append(['System', 'adjusting surface interval number'])
    args.extend(one_of_valid(config, ones))
    if (interval_adjustment == IntervalAdjustment.FOCAL_POSITION or
            interval_adjustment == IntervalAdjustment.FOCAL_POSITION_LENGTH or
            interval_adjustment ==
            IntervalAdjustment.FOCAL_POSITION_SEPARATED):
        args.append(['System', 'focal position'])
    if (interval_adjustment == IntervalAdjustment.TELEPHOTO_POWER_0):
        args.append(['Telephoto', 'system position'])
    if (interval_adjustment == IntervalAdjustment.FOCAL_LENGTH or
            interval_adjustment == IntervalAdjustment.FOCAL_POSITION_LENGTH):
        args.append(['System', 'focal length'])
    if (interval_adjustment == IntervalAdjustment.FOCAL_POSITION_SEPARATED):
        args.append(['System', 'fixed lens position'])
    if (interval_adjustment == IntervalAdjustment.FOCAL_POSITION_LENGTH):
        args.append(['System', 'iris position'])
    if (interval_adjustment == IntervalAdjustment.FOCAL_POSITION_SEPARATED):
        args.append(['System', 'adjusting interval number behind'])
    if (interval_adjustment == IntervalAdjustment.TELEPHOTO_POWER_0):
        ones = []
        ones.append(['Telephoto', 'number of surfaces'])
        ones.append(['Telephoto', 'number of powers'])
        args.extend(one_of_valid(config, ones))
    write_dataline(fout, config, args)

    fout.write('*NSPR:data ')
    args = []
    args.append(['System', 'surfaces of powers'])
    args.append(['System', 'division ratios of powers'])
    write_dataline(fout, config, args)

    fout.write('*INC:data ')
    args = []
    args.append([str(object_distance.number)])
    args.append([str(incident_light_range.number)])
    ones = []
    ones.append(['System', 'aperture radius'])
    ones.append(['System', 'light source radius'])
    args.extend(one_of_valid(config, ones))
    args.append(['System', 'object distance'])
    args.append([str(tracing_direction.number)])
    if (tracing_direction == TracingDirection.SCREEN_TO_APERTURE_EXT):
        args.append(['System', 'aperture to light source'])
        args.extend(all_of_valid(config, [['System', 'mirror to valve']]))
    elif (tracing_direction == TracingDirection.SOURCE_TO_APERTURE or
          tracing_direction == TracingDirection.SOURCE_TO_APERTURE_EXT):
        ones = []
        ones.append(['System', 'aperture position'])
        #        ones.append(['System', 'aperture radius real'])
        args.extend(all_of_valid(config, ones))
    write_dataline(fout, config, args)

    fout.write('*INCA:data ')
    args = []
    if (incident_light_range == IncidentLightRange.SYSTEM_SOURCE or
            incident_light_range == IncidentLightRange.APERTURE_IRIS_SOURCE):
        args.append(['Light Source', 'light source file'])
        args.append(['Light Source', 'diffusion radius'])
        args.append(['Light Source', 'number of filaments'])
        args.append(['Light Source', 'filament position'])
        args.append(['Light Source', 'filament length'])
        args.append(['Light Source', 'number of aperture segments'])
        args.append(['Light Source', 'angle range'])
        args.append(['Light Source', 'number of tracing rays'])
        args.append(['Light Source', 'iris surface number'])
        args.append(['Light Source', 'limit surface number'])
    write_dataline(fout, config, args)

    fout.write('*SCRN:data ')
    args = []
    if (tracing_direction == TracingDirection.SOURCE_TO_APERTURE):
        args.append(['System', 'light source to screen'])
        args.append(['System', 'screen radius'])
    write_dataline(fout, config, args)

    fout.write('*EF:data ')
    write_dataline(fout, config, [[
        str(input_data_type.number +
            (tracing_number.number + interval_adjustment.number) * 10)
    ]])

    fout.write('*DSCV:data ')
    args = []
    if (input_data_type == InputDataType.SURFACE):
        args.append(['System', 'surface data'])
    write_dataline(fout, config, args)

    fout.write('*DSPW:data ')
    args = []
    if (input_data_type == InputDataType.POWER):
        args.append(['System', 'power data'])
    elif (input_data_type == InputDataType.POWER_LIGHT):
        args.append(['System', 'power and light data'])
    elif (input_data_type == InputDataType.GROUP):
        args.append(['System', 'group data'])
    write_dataline(fout, config, args)

    fout.write('*PSPL:data ')
    args = []
    if (input_data_type == InputDataType.GROUP):
        for e in split_patterns:
            args.append(str(e.number))
        args.append(['System', 'split parameters first'])
        args.append(['System', 'split parameters second'])
    write_dataline(fout, config, args)

    fout.write('*GLAS:data ')
    write_dataline(fout, config, [['System', 'glass type names']])

    fout.write('*DSTK:data ')
    write_dataline(fout, config, [['System', 'surface intervals of powers']])

    fout.write('*CURV:data ')
    args = []
    for e in ptlc_specifications:
        args.append([str(e.number)])
    args.append(['System', 'power to lens curvature values'])
    write_dataline(fout, config, args)

    fout.write('*LRAD:data ')
    write_dataline(fout, config, [['System', 'lens radiuses']])

    fout.write('*LAR:data ')
    write_dataline(fout, config, [['System', 'lens effective radiuses']])

    fout.write('*ALR:data ')
    write_dataline(fout, config, [['System', 'surfaces axis length ratios']])

    fout.write('*ASPH:data ')
    write_dataline(fout, config, [['System',
                                   'surfaces aspheric coefficients']])

    fout.write('*DISP:data ')
    args = []
    args.append(['View Parameters', 'X-axis length'])
    args.append(['View Parameters', 'screen edge to first surface'])
    args.append(['View Parameters', 'Y-coordinate at screen center'])
    args.append(['View Parameters', 'X-axis grid pitch'])
    args.append(['View Parameters', 'Y-axis grid pitch'])
    args.append([str(screen_display_dimension.number)])
    write_dataline(fout, config, args)

    fout.write('*DEND:data ')
    args = []
    args.append(['Light Source', 'number of mirror splitting'])
    args.append(['Light Source', 'number of mirror facets'])
    args.append(['Light Source', 'heights of mirror facet edges'])
    args.append(
        ['Light Source', 'angles between mirror facets and optical axis'])
    write_dataline(fout, config, args)
