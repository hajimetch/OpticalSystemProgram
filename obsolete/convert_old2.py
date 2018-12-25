# -*- coding: utf-8 -*-
"""
Created on Fri Mar 16 09:59:09 2018

@author: hajimetch
"""

import configparser
from enum import Enum

#
# Associate setting items with each string expression.
#


class Projection(Enum):
    SPOT = 'spot'
    PATTERN = 'pattern'


class Reflector(Enum):
    TRUE = 'true'
    FALSE = 'false'


class Telephoto(Enum):
    TRUE = 'true'
    FALSE = 'false'


class InputDataType(Enum):
    SURFACE = 'surface data'
    POWER = 'power data'
    POWER_LIGHT = 'power and light data'
    GROUP = 'group data'


class TracingNumber(Enum):
    SINGLE = 'single'
    MULTI = 'multi'


class IntervalAdjustment(Enum):
    NONE = 'none'
    FOCAL_LENGTH = 'specify focal length'
    FOCAL_POSITION = 'specify focal position'
    TELEPHOTO_POWER_0 = 'telephoto system power set to 0'
    FOCAL_POSITION_LENGTH = 'specify focal position and length'


class ObjectDistance(Enum):
    INFINITY = 'to infinity'
    FIRST_POWER = 'to first power'
    FIRST_SURFACE = 'to first surface'
    PRINCIPAL_POSITION = 'to principal position of projection side'


class IncidentLightRange(Enum):
    SYSTEM_SOURCE = 'by system and light source'
    APERTURE_IRIS_SOURCE = 'by lens aperture iris and light source'
    APERTURE_IRIS = 'by lens aperture and iris'


class TracingDirection(Enum):
    SCREEN_TO_APERTURE = 'screen to aperture'
    SCREEN_TO_APERTURE_EXT = 'screen to aperture extended'
    SOURCE_TO_APERTURE = 'light source to aperture'
    SOURCE_TO_APERTURE_EXT = 'light source to aperture extended'


class SplitPattern(Enum):
    NONE = 'none'
    FIRST_INTERVAL = 'first power and interval'
    LAST_INTERVAL = 'last power and interval'
    RATIO_INTERVAL = 'power ratio and interval'
    FIRST_LAST = 'first and last powers'


class PTLCSpecification(Enum):
    FRONT = 'front'
    BACK = 'back'
    FRONT_OVER_BACK = 'front/back'
    BACK_OVER_FRONT = 'back/front'
    FIRST_LENS = 'first lens'
    SECOND_LENS = 'second lens'
    APERTURE_MIRROR = 'aperture or mirror'
    REFLECTOR = 'reflector'


class ScreenDisplayDimension(Enum):
    XYPLANE_COORDINATES = 'XY-plane and ray passing point coordinates'
    XYPLANE_3D = 'XY-plane and 3D'
    XYPLANE_XZPLANE = 'XY-plane and XZ-plane'
    SPLIT = 'split view'


#
# Number each setting item.
#

NUMBER = {}

NUMBER[Telephoto.TRUE] = 1
NUMBER[Telephoto.FALSE] = 0

NUMBER[InputDataType.SURFACE] = 0
NUMBER[InputDataType.POWER] = 1
NUMBER[InputDataType.POWER_LIGHT] = 2
NUMBER[InputDataType.GROUP] = 3

NUMBER[TracingNumber.SINGLE] = 0
NUMBER[TracingNumber.MULTI] = 5

NUMBER[IntervalAdjustment.NONE] = 0
NUMBER[IntervalAdjustment.FOCAL_LENGTH] = 1
NUMBER[IntervalAdjustment.FOCAL_POSITION] = 2
NUMBER[IntervalAdjustment.TELEPHOTO_POWER_0] = 3
NUMBER[IntervalAdjustment.FOCAL_POSITION_LENGTH] = 4

NUMBER[ObjectDistance.INFINITY] = 0
NUMBER[ObjectDistance.FIRST_POWER] = 1
NUMBER[ObjectDistance.FIRST_SURFACE] = 2
NUMBER[ObjectDistance.PRINCIPAL_POSITION] = 3

NUMBER[IncidentLightRange.SYSTEM_SOURCE] = 1
NUMBER[IncidentLightRange.APERTURE_IRIS_SOURCE] = 2
NUMBER[IncidentLightRange.APERTURE_IRIS] = 3

NUMBER[TracingDirection.SCREEN_TO_APERTURE] = 0
NUMBER[TracingDirection.SCREEN_TO_APERTURE_EXT] = -1
NUMBER[TracingDirection.SOURCE_TO_APERTURE] = 1
NUMBER[TracingDirection.SOURCE_TO_APERTURE_EXT] = 2

NUMBER[SplitPattern.NONE] = 0
NUMBER[SplitPattern.FIRST_INTERVAL] = 1
NUMBER[SplitPattern.LAST_INTERVAL] = 2
NUMBER[SplitPattern.RATIO_INTERVAL] = 3
NUMBER[SplitPattern.FIRST_LAST] = 4

NUMBER[PTLCSpecification.FRONT] = 0
NUMBER[PTLCSpecification.BACK] = 1
NUMBER[PTLCSpecification.FRONT_OVER_BACK] = 2
NUMBER[PTLCSpecification.BACK_OVER_FRONT] = 3
NUMBER[PTLCSpecification.FIRST_LENS] = 4
NUMBER[PTLCSpecification.SECOND_LENS] = 5
NUMBER[PTLCSpecification.APERTURE_MIRROR] = 6
NUMBER[PTLCSpecification.REFLECTOR] = 7

NUMBER[ScreenDisplayDimension.XYPLANE_COORDINATES] = 0
NUMBER[ScreenDisplayDimension.XYPLANE_3D] = 1
NUMBER[ScreenDisplayDimension.XYPLANE_XZPLANE] = 2
NUMBER[ScreenDisplayDimension.SPLIT] = 3

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
        raise ElementBoundsException(Exception)


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
config.read('examplemac.ini')

projection = Projection(
    config.get('Parameters', 'projection', fallback='spot'))
reflector = Reflector(config.get('Parameters', 'reflector', fallback='false'))
telephoto = Telephoto(config.get('Parameters', 'telephoto', fallback='false'))
input_data_type = InputDataType(
    config.get('Parameters', 'input data type', fallback='surface data'))
tracing_number = TracingNumber(
    config.get('Parameters', 'tracing number', fallback='single'))
interval_adjustment = IntervalAdjustment(
    config.get(
        'Parameters', 'lenses/powers interval adjustment', fallback='none'))
object_distance = ObjectDistance(
    config.get('Parameters', 'object distance', fallback='to infinity'))
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
        fallback='front').split(','):
    ptlc_specifications.append(PTLCSpecification(e.strip()))

with open('datafile.txt', 'wt') as fout:
    fout.write('*SDT:data ')
    sysname = 'FS' if projection == Projection.SPOT else ''
    sysname += 'REF' if reflector == Reflector.TRUE else ''
    sysname += config.get('Parameters', 'system name').strip()
    write_dataline(fout, config, [[sysname]])

    fout.write('*IAF:data ')
    write_dataline(fout, config, [[str(NUMBER[telephoto])]])

    fout.write('*MSDT:data ')
    args = []
    args.append(['Telephoto', 'distance between systems'])
    args.append(['Telephoto', 'focal length'])
    args.append(['Telephoto', 'principal position of object side'])
    args.append(['Telephoto', 'principal position of projection side'])
    args.append(['Telephoto', 'aperture for display'])
    args.append(['Telephoto', 'refractive index of projection region'])
    args.append(['Telephoto', 'adjusting lenses/powers interval number'])
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
    ones.append(['System', 'adjusting surface interval number'])
    ones.append(['System', 'adjusting power interval number'])
    args.extend(one_of_valid(config, ones))
    if (interval_adjustment == IntervalAdjustment.FOCAL_POSITION or
            interval_adjustment == IntervalAdjustment.FOCAL_POSITION_LENGTH):
        args.append(['System', 'focal position'])
    if (interval_adjustment == IntervalAdjustment.TELEPHOTO_POWER_0):
        args.append(['Telephoto', 'system position'])
    if (interval_adjustment == IntervalAdjustment.FOCAL_LENGTH or
            interval_adjustment == IntervalAdjustment.FOCAL_POSITION_LENGTH):
        args.append(['System', 'focal length'])
    if (interval_adjustment == IntervalAdjustment.FOCAL_POSITION_LENGTH):
        args.append(['System', 'iris position'])
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
    args.append([str(NUMBER[object_distance])])
    args.append([str(NUMBER[incident_light_range])])
    ones = []
    ones.append(['System', 'aperture radius'])
    ones.append(['System', 'light source radius'])
    args.extend(one_of_valid(config, ones))
    args.append(['System', 'object distance'])
    args.append([str(NUMBER[tracing_direction])])
    if (tracing_direction == TracingDirection.SCREEN_TO_APERTURE_EXT):
        args.append(['System', 'aperture to light source'])
        args.extend(all_of_valid(config, [['System', 'mirror to valve']]))
    elif (tracing_direction == TracingDirection.SOURCE_TO_APERTURE or
          tracing_direction == TracingDirection.SOURCE_TO_APERTURE_EXT):
        ones = []
        ones.append(['System', 'aperture position'])
        ones.append(['System', 'aperture radius real'])
        args.extend(all_of_valid(config, ones))
    write_dataline(fout, config, args)

    fout.write('*INCA:data ')
    args = []
    if (incident_light_range == IncidentLightRange.SYSTEM_SOURCE or
            incident_light_range == IncidentLightRange.APERTURE_IRIS_SOURCE):
        args.append(['Light Source', 'light source file'])
        args.append(['Light Source', 'diffusion radius'])
        args.append(['Light Source', 'filament number'])
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
        str(NUMBER[input_data_type] +
            (NUMBER[tracing_number] + NUMBER[interval_adjustment]) * 10)
    ]])

    fout.write('*DSCV:data ')
    args = []
    if (input_data_type == InputDataType.SURFACE):
        args.append(['System', 'surface intervals'])
        args.append(['System', 'surface curvature radiuses'])
    write_dataline(fout, config, args)

    fout.write('*DSPW:data ')
    args = []
    if (input_data_type == InputDataType.POWER):
        args.append(['System', 'power intervals'])
        args.append(['System', 'power values'])
    elif (input_data_type == InputDataType.POWER_LIGHT):
        args = []
        args.append(['System', 'incident heights at powers start'])
        args.append(['System', 'power intervals'])
        args.append(['System', 'incident heights at powers'])
    elif (input_data_type == InputDataType.GROUP):
        args = []
        args.append(['System', 'incident heights at groups start'])
        args.append(['System', 'group intervals'])
        args.append(['System', 'incident heights at groups'])
    write_dataline(fout, config, args)

    fout.write('*PSPL:data ')
    args = []
    if (input_data_type == InputDataType.GROUP):
        for e in split_patterns:
            args.append(str(NUMBER[e]))
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
        args.append([str(NUMBER[e])])
    args.append(['System', 'curvature specification values'])
    write_dataline(fout, config, args)

    fout.write('*LRAD:data ')
    write_dataline(fout, config, [['System', 'lens radiuses']])

    fout.write('*LAR:data ')
    write_dataline(fout, config, [['System', 'lens effective radiuses']])

    fout.write('*ALR:data ')
    write_dataline(fout, config, [['System', 'surfaces axis length ratios']])

    fout.write('*ASPH:data ')
    write_dataline(fout, config, [['System', 'surfaces aspheric coefficient']])

    fout.write('*DISP:data ')
    args = []
    args.append(['View Parameters', 'X-axis length'])
    args.append(['View Parameters', 'screen edge to first surface'])
    args.append(['View Parameters', 'Y-coordinate at screen center'])
    args.append(['View Parameters', 'X-axis grid pitch'])
    args.append(['View Parameters', 'Y-axis grid pitch'])
    args.append([str(NUMBER[screen_display_dimension])])
    write_dataline(fout, config, args)

    fout.write('*DEND:data ')
    args = []
    args.append(['Light Source', 'number of mirror splitting'])
    args.append(['Light Source', 'number of mirror facets'])
    args.append(['Light Source', 'heights of mirror facet edges'])
    args.append(
        ['Light Source', 'angles between mirror facets and optical axis'])
    write_dataline(fout, config, args)
