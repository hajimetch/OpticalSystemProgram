# -*- coding: utf-8 -*-
"""
Created on Thu Apr 05 10:22:12 2018

@author: hajimetch
"""

import configparser
import math
import numpy as np
from enum import Enum
from datetime import datetime

#
# Associate setting items with each string expression and number.
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
    FOCAL_POSITION_SEPARATED = 'specify focal position of separated system'


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
    NONE = 'none'
    FRONT = 'front'
    BACK = 'back'
    FRONT_OVER_BACK = 'front/back'
    BACK_OVER_FRONT = 'back/front'


class ScreenDisplayDimension(Enum):
    XYPLANE_COORDINATES = 'XY-plane and ray passing point coordinates'
    XYPLANE_3D = 'XY-plane and 3D'
    XYPLANE_XZPLANE = 'XY-plane and XZ-plane'
    SPLIT = 'split view'


class KindOfLens(Enum):
    SINGLE_LENS = 'single lens'
    COMBINED_LENS_FRONT = 'combined lens front'
    COMBINED_LENS_BACK = 'combined lens back'
    IRIS = 'iris'
    LS_REFLECTOR = 'light source reflector'
    REFLECTOR = 'reflector'
    APERTURE = 'aperture'
    LIGHT_SOURCE = 'light source'


class Direction(Enum):
    FORWARD = 1
    BACKWARD = -1


class OpticalSystem():
    def __init__(self, config):
        self.start_dt = datetime.now()

        self.read_data(config)
        self.check_data()

        self.prepare_lens_eff_radius(config)
        self.prepare_ptlc(config)
        self.prepare_axis_and_aspheric(config)
        self.prepare_glass_and_ri(config)

        if self.input_dtype == InputDataType.SURFACE:
            self.process_surf(config)
        elif self.input_dtype == InputDataType.POWER:
            self.process_power(config)
        elif self.input_dtype == InputDataType.POWER_LIGHT:
            self.process_power_light(config)
        elif self.input_dtype == InputDataType.GROUP:
            self.process_group(config)

        self.determine_hole_and_max()
        self.determine_kind_of_lens()

        self.calculate_inverse()

        if self.input_dtype != InputDataType.SURFACE and (
                self.interval_adj == IntervalAdjustment.FOCAL_LENGTH or
                self.interval_adj == IntervalAdjustment.TELEPHOTO_POWER_0):
            self.modify_lens_position()

        if (self.input_dtype == InputDataType.SURFACE and
                self.interval_adj == IntervalAdjustment.NONE):
            self.convert_surf_to_power()
        else:
            self.convert_power_to_surf()

    def read_data(self, config):
        """
        configで指定されたデータファイルから、光学系のデータを読み込む
        対応するBASICコード：*IN1 など
        """
        self.projection = Projection(
            config.get('Parameters', 'projection', fallback='spot'))
        self.telephoto = Telephoto(
            config.get('System', 'telephoto', fallback='false'))
        self.reflector = Reflector(
            config.get('System', 'reflector', fallback='false'))
        self.scale = config.getfloat('System', 'scale', fallback=100)
        self.input_dtype = InputDataType(
            config.get('System', 'input data type', fallback='surface data'))
        self.tracing_no = TracingNumber(
            config.get('System', 'tracing number', fallback='single'))
        self.interval_adj = IntervalAdjustment(
            config.get('System', 'interval adjustment', fallback='none'))
        self.obj_dist_spec = ObjectDistance(
            config.get(
                'Parameters',
                'object distance specification',
                fallback='to infinity'))
        self.incident_light_range = IncidentLightRange(
            config.get(
                'Parameters',
                'incident light range definition',
                fallback='by system and light source'))
        self.tracing_dir = TracingDirection(
            config.get(
                'Parameters',
                'tracing direction',
                fallback='screen to aperture'))

        self.system_name = config.get('System', 'system name', fallback='')
        self.no_powers = config.getint(
            'System', 'number of powers', fallback=0)
        self.no_surfs = config.getint(
            'System', 'number of surfaces', fallback=0)
        self.no_lenses = config.getint(
            'System', 'number of lenses', fallback=0)
        self.no_groups = config.getint(
            'System', 'number of groups', fallback=self.no_powers)
        self.adj_power_inr_no = config.getint(
            'System', 'adjusting power interval number', fallback=0)
        self.adj_surf_inr_no = config.getint(
            'System', 'adjusting surface interval number', fallback=0)
        self.adj_inr_no_bhnd = config.getint(
            'System', 'adjusting interval number behind', fallback=0)
        self.focal_pos = config.getfloat(
            'System', 'focal position', fallback=0) / self.scale
        self.focal_len = config.getfloat(
            'System', 'focal length', fallback=0) / self.scale
        self.iris_pos = config.getfloat(
            'System', 'iris position', fallback=0) / self.scale
        self.fixed_lens_pos = config.getfloat(
            'System', 'fixed lens position', fallback=0) / self.scale
        self.ap_radius = config.getfloat(
            'System', 'aperture radius', fallback=0) / self.scale
        self.ls_radius = config.getfloat(
            'System', 'light source radius', fallback=0) / self.scale
        self.obj_dist = config.getfloat(
            'System', 'object distance', fallback=0) / self.scale
        self.ap_to_ls = config.getfloat(
            'System', 'aperture to light source', fallback=0) / self.scale
        self.mirror_to_valve = config.getfloat(
            'System', 'mirror to valve', fallback=0) / self.scale

        self.tel_no_powers = config.getint(
            'Telephoto', 'number of powers', fallback=0)
        self.tel_no_surfs = config.getint(
            'Telephoto', 'number of surfaces', fallback=0)
        self.tel_system_pos = config.getfloat(
            'Telephoto', 'system position', fallback=0) / self.scale
        self.tel_system_dist = config.getfloat(
            'Telephoto', 'distance between systems', fallback=0) / self.scale
        self.tel_focal_len = config.getfloat(
            'Telephoto', 'focal length', fallback=0) / self.scale
        self.tel_prin_pos_obj = config.getfloat(
            'Telephoto', 'principal position of object side',
            fallback=0) / self.scale + self.tel_system_dist
        self.tel_prin_pos_proj = config.getfloat(
            'Telephoto', 'principal position of projection side',
            fallback=0) / self.scale + self.tel_system_dist
        self.tel_ap_disp = config.getfloat(
            'Telephoto', 'aperture for display', fallback=0) / self.scale
        self.tel_ri_proj = config.getfloat(
            'Telephoto', 'refractive index of projection region', fallback=0)
        self.tel_adj_inr_no = config.getint(
            'Telephoto', 'adjusting interval number', fallback=0)

        self.ls_file = config.get(
            'Light Source', 'light source file', fallback='')
        self.ls_diff_radius = config.getfloat(
            'Light Source', 'diffusion radius', fallback=0) / self.scale
        self.ls_no_filamets = config.getint(
            'Light Source', 'number of filaments', fallback=0)
        self.ls_filament_pos = config.getfloat(
            'Light Source', 'filament position', fallback=0) / self.scale
        self.ls_filament_len = config.getint(
            'Light Source', 'filament length', fallback=0) / self.scale
        self.ls_no_segments = config.getint(
            'Light Source', 'number of aperture segments', fallback=0)
        self.ls_angle_range = config.getfloat(
            'Light Source', 'angle range', fallback=0) / self.scale
        self.ls_no_rays = config.getint(
            'Light Source', 'number of tracing rays', fallback=0)
        self.ls_iris_no = config.getint(
            'Light Source', 'iris surface number', fallback=0)
        self.ls_limit_no = config.getint(
            'Light Source', 'limit surface number', fallback=0)

        self.surfs_of_power = get_intlist(config, 'System',
                                          'surfaces of powers')
        self.lens_radius = get_floatlist(config, 'System', 'lens radius')
        self.projection_dist = get_floatlist(config, 'System',
                                             'projection distances')
        self.surf_inr_of_power = get_floatlist(config, 'System',
                                               'surface intervals of powers')
        self.division_ratio = get_floatlist(config, 'System',
                                            'division ratios of powers')
        self.split_param_1 = get_floatlist(config, 'System',
                                           'split parameters first')
        self.split_param_2 = get_floatlist(config, 'System',
                                           'split parameters second')

    def check_data(self):
        """
        read_data()で読み込んだデータに不整合がないかチェックする
        """
        if (self.telephoto == Telephoto.FALSE and
                self.interval_adj == IntervalAdjustment.TELEPHOTO_POWER_0):
            raise DataConsistenceError(
                '[System: telephoto] is "false"' +
                ' though [System: interval adjustment] is' +
                ' "telephoto system power set to 0".')
        if (self.surfs_of_power is not None and
                self.input_dtype == InputDataType.SURFACE and
                self.interval_adj == IntervalAdjustment.NONE):
            raise DataConsistenceError(
                '[System: surface intervals of powers] is not void' +
                ' though [System: input data type] is' +
                ' "surface" and [System: interval adjustment] is "none"')

    def prepare_lens_eff_radius(self, config):
        """
        configで指定されたデータファイルのlens effective radiusの項目を展開して
        lens_eff_radiusリストを構成する
        対応するBASICコード：*IN6
        """
        lens_eff_radius_data = get_floatlist(config, 'System',
                                             'lens effective radius')
        lens_eff_radius_dict = {
            i: {
                'front': j,
                'back': k
            }
            for i, j, k in zip(* [iter(lens_eff_radius_data)] * 3)
        }
        self.lens_eff_radius = []
        for lens_no in range(0, self.no_lenses):
            if lens_no not in lens_eff_radius_dict:
                self.lens_eff_radius.append(self.lens_radius[lens_no])
                self.lens_eff_radius.append(self.lens_radius[lens_no])
            else:
                if (0 < lens_eff_radius_dict[lens_no]['front'] and
                        lens_eff_radius_dict[lens_no]['front'] <
                        self.lens_radius[lens_no]):
                    self.lens_eff_radius.append(
                        lens_eff_radius_dict[lens_no]['front'])
                else:
                    self.lens_eff_radius.append(self.lens_radius[lens_no])
                if (0 < lens_eff_radius_dict[lens_no]['back'] and
                        lens_eff_radius_dict[lens_no]['back'] <
                        self.lens_radius[lens_no]):
                    self.lens_eff_radius.append(
                        lens_eff_radius_dict[lens_no]['back'])
                else:
                    self.lens_eff_radius.append(self.lens_radius[lens_no])

    def prepare_ptlc(self, config):
        """
        configで指定されたデータファイルのpower to lens curvature specificationsの項目を展開して
        ptlc_specリストとptlc_valリストを構成する
        対応するBASICコード：*IN6
        """
        self.ptlc_spec = []
        for e in config.get(
                'System',
                'power to lens curvature specifications',
                fallback='none').split(','):
            self.ptlc_spec.append(PTLCSpecification(e.strip()))
            self.ptlc_val = config.get_floatlist(config, 'System',
                                                 'power to lens values')

    def prepare_axis_and_aspheric(self, config):
        """
        configで指定されたデータファイルのsurfaces axis length ratiosと
        surfaces aspheric coefficientsの項目を読み込み、
        axis_ratioリストとaspheric_coefficientリストを構成する
        対応するBASICコード：*IN6
        """
        surf_axis_data = get_floatlist(config, 'System',
                                       'surface axis length ratios')
        surf_axis_dict = dict(zip(* [iter(surf_axis_data)] * 2))
        surf_aspheric_data = get_floatlist(config, 'System',
                                           'surfaces aspheric coefficients')
        self.axis_ratio = []
        self.aspheric_coefficient = []
        aspheric_no = 0
        for surf_no in range(0, self.no_surfs):
            if surf_no not in surf_axis_dict:
                self.axis_ratio.append(1)
                self.aspheric_coefficient.append([])
            else:
                self.axis_ratio.append(surf_axis_dict[surf_no])
                if surf_axis_dict[surf_no] >= 0:
                    self.aspheric_coefficient.append([])
                else:
                    aspheric_coefficient_element = []
                    for index in range(0, 4):
                        aspheric_coefficient_element.append(
                            surf_aspheric_data[aspheric_no * 4 + index] *
                            (self.scale / 100)**(index * 2 + 1))
                        self.aspheric_coefficient.append(
                            aspheric_coefficient_element)
                        aspheric_no += 1

    def prepare_glass_and_ri(self, config):
        """
        予め用意されたGLASS_DATAおよび
        configで指定されたデータファイルのglass type namesの項目をもとに
        光学系を定める以下の各リストを構成する
        glass, direction, dispersion, ri, signed_ri
        対応するBASICコード：*IN6
        """
        global GLASS_DATA
        dr = 1
        self.glass = []
        self.direction = []
        self.dispersion = []
        self.ri = []
        self.signed_ri = []
        for glass_name in get_stringlist(config, 'System', 'glass type names'):
            glass_name = 'AIR' if glass_name == '' else glass_name
            if glass_name == 'REF':
                dr = -dr
            if glass_name not in GLASS_DATA:
                raise NoGlassDataError(glass_name)
            self.glass.append(GLASS_DATA[glass_name])
            self.direction.append(Direction(dr))
            self.dispersion.append(GLASS_DATA[glass_name]['dispersion'])
            self.ri.append(GLASS_DATA[glass_name]['nd'])
            self.signed_ri.append(dr * GLASS_DATA[glass_name]['nd'])
            self.diff_ri = [
                j - i for i, j in zip(* [iter(self.signed_ri)] * 2)
            ]

    def process_surf(self, config):
        """
        configで指定されたデータファイルのsurface dataの項目から
        surf_inrリストとsurf_rocリストを読み込む
        対応するBASICコード：*IN6
        """
        surf_data = get_floatlist(config, 'System', 'surface data')
        self.surf_inr = surf_data[::2]
        self.surf_roc = surf_data[1::2]

    def process_power(self, config):
        """
        configで指定されたデータファイルのpower dataの項目から
        power_inrリストとpower_valリストを読み込む
        対応するBASICコード：*IN6
        """
        power_data = get_floatlist(config, 'System', 'power data')
        self.power_inr = power_data[::2]
        self.power_val = power_data[1::2]
        self.power_inr_0 = self.power_inr[0]

    def process_power_light(self, config):
        """
        configで指定されたデータファイルのpower and light dataの項目から
        light_hgt_at_powerリストとpower_inrリストを読み込み、
        power_valリストの値を計算する
        対応するBASICコード：*IN6, *EF2
        """
        power_light_data = get_floatlist(config, 'System',
                                         'power and light data')
        self.light_hgt_at_power = power_light_data[::2]
        self.power_inr = power_light_data[1::2]
        surf_no = 0
        self.power_val = []
        for power_no in range(0, self.no_powers):
            if (self.diff_ri[surf_no] == 0 and
                    self.surfs_of_power[power_no] == 1 or
                    self.power_inr[power_no] == 0):
                self.power_val.append(0)
            else:
                self.power_val.append(
                    ((self.light_hgt_at_power[power_no + 1] - self.
                      light_hgt_at_power[power_no + 2]
                      ) * rev(self.power_inr[power_no + 1]) -
                     (self.light_hgt_at_power[power_no] - self.
                      light_hgt_at_power[power_no +
                                         1]) * rev(self.power_inr[power_no])) /
                    self.light_hgt_at_power[power_no + 1])
                surf_no += self.surfs_of_power[power_no]
        self.power_inr_0 = self.group_inr[0]

    def process_group(self, config):
        """
        configで指定されたデータファイルのgroup dataの項目から
        light_hgt_at_groupリストとgroup_inrリストを読み込み、
        split patternsの項目から
        split_patternリストを読み込む。
        group_valリストとpower_inrリストとpower_valリストを計算する
        対応するBASICコード：*IN6, *EF3
        """
        group_data = get_floatlist(config, 'System', 'group data')
        self.light_hgt_at_group = group_data[::2]
        self.group_inr = group_data[1::2]
        self.split_pattern = [
            SplitPattern(e.strip())
            for e in get_stringlist(config, 'System', 'split patterns')
        ]

        # light_hgt_at_groupおよびgroup_inrからgroup_valを計算
        surf_no = 0
        power_no = 0
        self.group_val = []
        for group_no in range(0, self.no_groups):
            if (self.split_pattern[group_no] == SplitPattern.NONE and
                    self.diff_ri[surf_no] == 0 and
                    self.surfs_of_power[power_no] == 1 or self.group_inr == 0):
                self.group_val.append(0)
            else:
                self.group_val.append(
                    ((self.light_hgt_at_group[group_no + 1] - self.
                      light_hgt_at_group[group_no + 2]
                      ) * rev(self.group_inr[group_no + 1]) -
                     (self.light_hgt_at_group[group_no] - self.
                      light_hgt_at_group[group_no +
                                         1]) * rev(self.group_inr[group_no])) /
                    self.light_hgt_at_group[group_no + 1])
                power_no += 1
                surf_no += self.surfs_of_power[power_no]
            if self.split_pattern[group_no] != SplitPattern.NONE:
                power_no += 1
                surf_no += self.surfs_of_power[power_no]

        # group_inrおよびgroup_valからpower_inrとpower_valへ変換
        power_no = 0
        delta = 0
        self.power_inr = []
        self.power_val = []
        for group_no in range(0, self.no_groups):
            if self.split_pattern[group_no] == SplitPattern.NONE:
                self.power_inr.append(self.group_inr[group_no] - delta)
                self.power_val.append(self.group_val[group_no])
                delta = 0
                power_no += 1
            else:
                v_u = self.spliv_param_1[group_no]
                v_v = self.spliv_param_2[group_no]
                v_f = self.group_val[group_no]
                v_g = self.group_inr[group_no]
                if self.split_pattern[group_no] == SplitPattern.FIRST_INTERVAL:
                    self.power_inr.append(v_g - v_v * (v_f - v_u) /
                                          (1 - v_u * v_v) / v_f - delta)
                    self.power_val.append(v_u)
                    self.power_inr.append(v_v)
                    self.power_val.append((v_f - v_u) / (1 - v_u * v_v))
                elif self.split_pattern[
                        group_no] == SplitPattern.LAST_INTERVAL:
                    self.power_inr.append(v_g - (v_u * v_v) / v_f - delta)
                    self.power_val.append((v_f - v_u) / (1 - v_u * v_v))
                    self.power_inr.append(v_v)
                    self.power_val.append(v_u)
                elif self.split_pattern[
                        group_no] == SplitPattern.RATIO_INTERVAL:
                    self.power_inr.append(
                        v_g - v_v * v_u / v_f * fn_p(v_u, v_f, v_v) - delta)
                    self.power_val.append(fn_p(v_u, v_f, v_v))
                    self.power_inr.append(v_v)
                    self.power_val.append(v_u * fn_p(v_u, v_f, v_v))
                elif self.split_pattern[group_no] == SplitPattern.FIRST_LAST:
                    self.power_inr.append(v_g - (v_u + v_v - v_f) / v_v -
                                          delta)
                    self.power_val.append(v_u)
                    self.power_inr.append((v_u + v_v - v_f) / v_u / v_v)
                    self.power_val.append(v_v)
                    delta = self.power_val[power_no] * self.power_inr[
                        power_no + 1] / v_f
                power_no += 2
        self.power_inr_0 = self.power_inr[0]

    def determine_hole_and_max(self):
        """
        lens_eff_radiusからhole_valueを
        lens_radiusとap_radiusからmax_radiusを
        それぞれ決定する
        対応するBASICコード：*IN6
        """
        self.hole_value = self.lens_eff_radius[self.no_lenses * 2 - 1]**2 if (
            self.tracing_dir == TracingDirection.SCREEN_TO_APERTURE_EXT and
            self.lens_eff_radius[self.no_lenses * 2 - 1] <
            self.lens_radius[self.no_lenses - 1]) else 0
        self.max_radius = max(self.lens_radius + [self.ap_radius])

    def determine_kind_of_lens(self):
        """
        光学系の構成要素の種類を順に決定し、kind_of_lensリストに格納する
        対応するBASICコード：*DRV
        """
        self.lens_index = []
        self.surf_eff_radius = []
        self.kind_of_lens = []
        surf_no = 0
        lens_no = 0
        for power_no in range(0, self.no_powers):
            if (self.adj_surf_inr_no == 0 and
                    self.adj_power_inr_no == power_no):
                self.adj_surf_inr_no = surf_no
            if (self.adj_power_inr_no == 0 and
                    self.adj_surf_inr_no == surf_no):
                self.adj_power_inr_no = power_no
            for surf in range(surf_no,
                              surf_no + self.surfs_of_power[power_no]):
                self.lens_index.append(lens_no)
                if self.is_air(surf) and self.is_air(surf + 1):
                    if (self.direction[surf].value ==
                            -self.direction[surf + 1].value):
                        self.kind_of_lens.append(KindOfLens.REFLECTOR)
                    elif (lens_no == self.no_lenses - 1 and self.tracing_dir ==
                          TracingDirection.SCREEN_TO_APERTURE_EXT) or (
                              surf == 0 and self.surf_inr[0] < 0 and
                              (self.tracing_dir == TracingDirection.
                               SOURCE_TO_APERTURE or self.tracing_dir ==
                               TracingDirection.SOURCE_TO_APERTURE_EXT)):
                        self.kind_of_lens.append(KindOfLens.LS_REFLECTOR)
                    else:
                        self.kind_of_lens.append(KindOfLens.IRIS)
                    self.surf_roc[surf] = 0
                    self.axis_ratio[surf] = 1
                    self.surf_eff_radius.append(
                        self.lens_eff_radius[lens_no * 2])
                    lens_no += 1
                elif self.is_air(surf) and not self.is_air(surf + 1):
                    self.kind_of_lens.append(KindOfLens.SINGLE_LENS)
                    self.surf_eff_radius.append(
                        self.lens_eff_radius[lens_no * 2])
                    lens_no += 1
                elif not self.is_air(surf) and not self.is_air(surf + 1):
                    self.kind_of_lens[-1] = KindOfLens.COMBINED_LENS_FRONT
                    self.kind_of_lens.append(KindOfLens.COMBINED_LENS_BACK)
                    self.surf_eff_radius.append(
                        min(self.lens_eff_radius[lens_no * 2 - 1],
                            self.lens_eff_radius[lens_no * 2]))
                    lens_no += 1
                else:
                    self.surf_eff_radius.append(
                        self.lens_eff_radius[lens_no * 2 - 1])
            surf_no = surf_no + self.surfs_of_power[power_no]
        self.surf_inr_0 = self.surf_inr[0]
        if self.tracing_dir == TracingDirection.SCREEN_TO_APERTURE_EXT:
            self.kind_of_lens.append(KindOfLens.APERTURE)
        elif (self.tracing_dir == TracingDirection.SOURCE_TO_APERTURE or
              self.tracing_dir == TracingDirection.SOURCE_TO_APERTURE_EXT):
            self.kind_of_lens.append(KindOfLens.LIGHT_SOURCE)
            self.kind_of_lens.append(KindOfLens.APERTURE)

    def calculate_inverse(self):
        """
        収差の計算に用いる逆行列を算出する
        対応するBASICコード：*IN5
        """
        CAL_POINTS = [.4, .6, .8, 1, 1.05],
        CAL_POINTS_EQUAL_AREA = [.020, .316, .548, .707, .837, .949, 1]

        self.inverse = {}
        tmp = []
        for index in range(4, 17, 2):
            tmp.append(sum([e**index for e in CAL_POINTS]))
        self.inverse['even'] = np.linalg.inv(
            np.matrix([[tmp[i + j] for i in range(0, 3)]
                       for j in range(0, 3)]))
        tmp = []
        for index in range(2, 11):
            tmp.append(sum([e**index for e in CAL_POINTS[::-1] + CAL_POINTS]))
        self.inverse['odd'] = np.linalg.inv(
            np.matrix([[tmp[i + j] for i in range(0, 4)]
                       for j in range(0, 4)]))
        tmp = []
        for index in range(4, 17, 2):
            tmp.append(sum([e**index for e in CAL_POINTS_EQUAL_AREA]))
        self.inverse['equal_area'] = np.linalg.inv(
            np.matrix([[tmp[i + j] for i in range(0, 3)]
                       for j in range(0, 3)]))

    def modify_lens_position(self):
        """
        絞りを含むレンズ位置の修正を行う
        対応するBASICコード：*INM
        """
        if self.interval_adj == IntervalAdjustment.FOCAL_LENGTH:
            dist_infinity = True
            telephoto = False
            rev_fl = self.signed_ri[self.no_surfs] * rev(self.focal_len)
        else:
            dist_infinity = (self.obj_dist_spec == ObjectDistance.INFINITY)
            telephoto = True
            rev_fl = 0
        v_y, v_a, v_p = self.ray_tracing(dist_infinity, telephoto)
        y1 = v_a - rev_fl
        x2 = self.power_inr[self.adj_power_inr_no]
        x1 = x2 * 1.1
        self.power_inr[self.adj_power_inr_no] = x1
        while True:
            v_y, v_a, v_p = self.ray_tracing(dist_infinity, telephoto)
            y2 = v_a - rev_fl
            x1, y1, x2, y2 = newton_step(x1, y1, x2, y2)
            self.power_inr[self.adj_power_inr_no] = x1
            if abs(y2) < .00001:
                break
        if self.interval_adj == IntervalAdjustment.TELEPHOTO_POWER_0:
            self.power_inr.append(self.tel_system_pos - v_p)

    def convert_surf_to_power(self):
        """
        屈折面データから、パワーデータへの変換を行う
        対応するBASICコード：*EFG
        """
        surf_no = 0
        self.power_inr = []
        self.power_val = []
        power_inr_prep_val = 0
        for power_no in range(0, self.no_powers):
            v_a = 0
            v_y = 1
            v_w = 0
            v_e_lst = []
            v_f_lst = []
            for surf in range(surf_no,
                              surf_no + self.surfs_of_power[power_no]):
                self.surf_inr_of_power.append(self.surf_inr[surf])
                v_e = self.surf_inr[surf] / self.signed_ri[surf]
                v_w = v_w + v_e / (v_y - v_e * v_a) / v_y
                v_y = v_y - v_e * v_a
                v_f = self.diff_ri[surf] * rev(self.surf_roc[surf])
                v_a = v_a + v_y * v_f
                v_e_lst.append(v_e)
                v_f_lst.append(v_f)
                lens_no = self.lens_index[surf]
                if (self.direction[surf].value == -self.direction[surf +
                                                                  1].value):
                    self.ptlc_val[lens_no] = self.surf_roc[surf]
                elif not self.is_air(surf + 1):
                    if (self.ptlc_spec[lens_no] == PTLCSpecification.FRONT):
                        self.ptlc_val[lens_no] = self.surf_roc[surf]
                    elif (self.ptlc_spec[lens_no] == PTLCSpecification.BACK):
                        self.ptlc_val[lens_no] = self.surf_roc[surf + 1]
                    elif (self.ptlc_spec[lens_no] ==
                          PTLCSpecification.FRONT_OVER_BACK or
                          self.ptlc_spec[lens_no] ==
                          PTLCSpecification.BACK_OVER_FRONT):
                        if abs(self.surf_roc[surf + 1]) < abs(
                                self.surf_roc[surf]):
                            self.ptlc_spec[
                                lens_no] = PTLCSpecification.BACK_OVER_FRONT
                            self.ptlc_val[lens_no] = (
                                self.surf_roc[surf +
                                              1] * rev(self.surf_roc[surf]))
                        else:
                            self.ptlc_spec[
                                lens_no] = PTLCSpecification.FRONT_OVER_BACK
                            self.ptlc_val[lens_no] = (self.surf_roc[
                                surf] * rev(self.surf_roc[surf + 1]))
                elif self.is_air(surf):
                    self.ptlc_val[lens_no] = self.power_val.append(v_a)
            v_a = rev(v_a)
            self.power_inr.append(power_inr_prep_val + (1 - 1 / v_y) * v_a +
                                  v_w)
            power_inr_prep_val = (1 - v_y) * v_a
            if self.surfs_of_power[power_no] == 2 and self.is_air(surf_no + 1):
                self.division_ratio.append(v_f_lst[1] / v_f_lst[0])
            elif self.surfs_of_power[power_no] == 3:
                if self.is_air(surf_no + 1):
                    self.division_ratio.append(
                        (v_f_lst[1] + v_f_lst[2] - v_e_lst[2] * v_f_lst[1] *
                         v_f_lst[2]) * rev(v_f_lst[0]))
                elif self.is_air(surf_no + 2):
                    self.division_ratio.append(
                        v_f_lst[2] / (v_f_lst[0] + v_f_lst[1] -
                                      v_e_lst[1] * v_f_lst[0] * v_f_lst[1]))
                else:
                    self.ptlc_val[self.lens_index[surf_no + 2] -
                                  2] = self.surf_roc[surf_no + 1]
                    if (self.direction[surf_no +
                                       1].value == self.direction[surf_no +
                                                                  2].value):
                        self.division_ratio.append(
                            -((self.signed_ri[surf_no +
                                              2] - self.direction[surf_no + 1].
                               value) * rev(self.surf_roc[surf_no + 1]) *
                              (1 - v_e_lst[2] * v_f_lst[2]) + v_f_lst[2]) /
                            ((self.signed_ri[surf_no +
                                             1] - self.direction[surf_no + 1].
                              value) * rev(self.surf_roc[surf_no + 1]) *
                             (1 - v_e_lst[1] * v_f_lst[0]) + v_f_lst[0]))
                    else:
                        self.division_ratio.append(
                            -(self.signed_ri[surf_no + 2] * rev(
                                self.surf_roc[surf_no + 1]) *
                              (1 - v_e_lst[2] * v_f_lst[2]) + v_f_lst[3]) /
                            (self.signed_ri[surf_no + 1] * rev(
                                self.surf_roc[surf_no + 1]) *
                             (1 - v_e_lst[1] * v_f_lst[0]) + v_f_lst[0]))
            elif self.surfs_of_power[power_no] == 4:
                self.division_ratio.append(
                    (v_f_lst[2] + v_f_lst[3] -
                     v_e_lst[3] * v_f_lst[2] * v_f_lst[3]) /
                    (v_f_lst[0] + v_f_lst[1] -
                     v_e_lst[1] * v_f_lst[0] * v_f_lst[1]))
            else:
                self.division_ratio.append(0)
            surf_no += self.surfs_of_power[power_no]

        if (self.obj_dist_spec == ObjectDistance.FIRST_POWER or
                self.obj_dist_spec == ObjectDistance.INFINITY and
                self.input_dtype != InputDataType.SURFACE and
                self.interval_adj != IntervalAdjustment.NONE):
            self.surf_inr[0] = (self.power_inr_0 - self.power_inr[0]
                                ) * self.signed_ri[0] + self.surf_inr_0
            self.power_inr[0] = self.power_inr_0
        else:
            self.power_inr_0 = self.power_inr[0]

    def convert_power_to_surf(self):
        """
        パワーデータから、屈折面データへの変換を行う
        対応するBASICコード：*DRG
        """
        surf_no = 0
        cal_delta = False
        LU = []
        for power_no in range(0, self.no_powers):
            lens_no = self.lens_index[surf_no]
            S = self.ptlc_spec[lens_no]
            if (self.ptlc_val[lens_no] != 0 and
                    self.ptlc_spec[lens_no] == PTLCSpecification.FRONT or
                    self.ptlc_spec[lens_no] == PTLCSpecification.BACK):
                R = -1 / self.ptlc_val[lens_no]
            else:
                R = self.ptlc_val[lens_no]
            L = []
            ED = []
            for surf in range(0, self.surfs_of_power[power_no]):
                if (1 < surf and 1 < self.ri[surf_no + surf] and
                        self.surf_inr_of_power[surf_no + surf] <= 0):
                    cal_delta = True
                L.append(self.diff_ri[surf_no + surf])
                ED.append(
                    abs(self.surf_inr_of_power[surf_no + surf]) /
                    self.signed_ri[surf_no + surf])
                LU.append(abs(self.surf_inr_of_power[surf_no + surf]))
            ED[0] = 0
            L0 = L[0]
            L1 = L[1]
            N1 = self.ri[surf_no + 1]
            while True:
                E1 = ED[1]
                if self.surfs_of_power[power_no] > 2:
                    L2 = L[2]
                    L3 = L[3]
                    N2 = self.ri[surf_no + 2]
                    E2 = ED[2]
                    E3 = ED[3]
                if self.surfs_of_power[power_no] == 1:
                    F0 = self.power_val[power_no]
                elif self.surfs_of_power[power] == 2:
                    if self.ri[surf_no + 1] == 1:
                        F0 = fn_p(self.division_ratio[power_no],
                                  self.power_val[power_no], ED[1])
                        F1 = F0 * self.division_ratio[power_no]
                    else:
                        F0, F1 = fn_ps(S, L0, L1, self.power_val[power_no], E1,
                                       R, F0, F1)
                elif self.surfs_of_power[power] == 3:
                    if self.ri[surf_no + 1] == 1:
                        pass
        """
        self.surf_inr[0] = 0
        I:surf_no=1
        for IT:power_no=1 to self.no_powers:
            J=I:surf_no
            KS=self.surfs_of_power[power_no]
            A=self.division_ratio[power_no]
            P0=self.power_val[power_no]
            IL:lens_no=self.lens_index[surf_no]
            for M:surf=1 to self.surfs_of_power[power_no]:
                Q=self.surf_inr_of_power[surf_no]
                if 1<surf and 1<self.ri[surf_no] and self.surf_inr_of_power[surf_no]<=0:
                    ILT:cal_delta=True
                L[surf]=LL[surf_no]
                Q=abs(Q)
                LU[surf_no]=Q
                ED[surf]=Q/self.signed_ri[surf_no]
                surf_no=surf_no+1
            ED[1]=0
            L1=L[1]
            L2=L[2]
            N2=self.ri[J+1]
            S=self.ptlc_spec[IL:lens_no]
            R=self.ptlc_val[IL:lens_no]
            if R!=0 and S=PTLCSpecification.FRONT or S=PTLCSpecification.BACK:
                R=-1/R
            while True:
                E2=ED[2]
                if KS:self.surfs_of_power[power_no]>2:
                    L3=L[3]
                    L4=L[4]
                    N3=NA[J+2]
                    E3=ED[3]
                    E4=ED[4]
                if KS:self.surfs_of_power[power_no]=1:
                    F1=P0:self.power_val[power_no]
                elif KS:self.surfs_of_power[power_no]=2:
                    if N2:self.ri[J+1]=1:
                        F1=fnP(A:self.division_ratio[power_no],P0:self.power_val[power_no],E2:ED[2])
                        F2=F1*A:self.division_ratio[power_no]
                    else:
                        F1,F2=PS(S:ptlc_spec,L[1],L[2],P0:self.power_val,ED[2],R,,)
                elif KS:self.surfs_of_power[power_no]=3:
                    if N2:self.ri[J+1]=1:
                        gosub *FC
                        F1=P1
                    elif N3:NA[J+2]=1:
                        gosub *FC
                        F3=P2
                    else:
                        R=rev(self.ptlc_val[self.lens_index[surf_no]+1])
                        if self.signed_ri[J+1]*self.signed_ri[J+2])>0:
                            O=sgn(self.signed_ri[J+1])
                        gosub *FC
                        F2=L2*R
                elif KS:self.surfs_of_power[power_no]=4:
                    S2=self.ptlc_spec[self.lens_index[surf_no]+1]
                    R2=self.ptlc_val[self.lens_index[surf_no]+1]
                    if R2!=0 and S2=PTLCSpecification.FRONT or S2=PTLCSpecification.Back:
                        R2=-1/R2
                    if N3:NS[J+2]=1:
                        gosub *FC
                    else:
                        ERROR
                else:
                    ERROR
                FD[1]=F1
                FD[2]=F2
                FD[3]=F3
                FD[4]=F4
                I:surf_no=J
                if ILT:cal_delta=False:
                    break
                else:
                    ILT=False
                Q=abs(self.surf_inr_of_power[surf_no])
                for M:surf=1 to KS:self.surfs_of_power[power_no]
                    B=Q
                    B=QR
                    F=FD[surf]
                    IL=ILD:self.lens_index[surf_no] - 1
                    Q=0
                    QR=0
                    if abs(F)>.000001 and LL[surf_no] != 0:
                        Q=LL:self.diff_ri[surf_no]/F/LC:self.axis_ratio[surf_no]**2
                        QR=1/R
                    if surf>0 and LT:self.surf_inr_of_power[surf_no]=<0 and 1<NA:self.ri[surf_no]:
                        H=LHD:self.lens_radius[IL]
                        O=-LT:self.surf_inr_of_power[surf_no]
                        if
        """

    def ray_tracing(self, dist_infinity, telephoto):
        v_y = 0 if dist_infinity else 1
        v_a = 1 - v_y
        v_p = -self.power_inr[0]
        for index in range(0, self.tel_no_powers
                           if telephoto else self.no_powers):
            v_y = v_y - v_a * self.power_inr[index]
            v_a = v_a + v_y * self.power_val[index]
            v_p = v_p + self.power_inr[index]
        return v_y, v_a, v_p

    def is_air(self, surf_no):
        """
        光学系内のsurf_noで示される場所が
        空気層か否かを返す関数
        """
        return (self.ri[surf_no] == 1 or
                (surf_no not in range(0, self.no_surfs)))


class DataConsistenceError(Exception):
    pass


class NoGlassDataError(Exception):
    pass


class IllegalNumberError(Exception):
    pass


def fn_fc(dri_lst, ri_lst, pv_lst, pi_lst, ptlc_v_lst, ptlc_s_lst, pv, divr,
          delta, no_surf, sign):
    """
    対応するBASICコード：*FC
    """
    d0 = 0.05
    d1 = 0
    t_pv_lst, pv_lst, q0 = fn_pc(dri_lst, ri_lst, pv_lst, pi_lst, ptlc_v_lst,
                                 ptlc_s_lst, pv, divr, d1, no_surf, sign)
    t_pv_lst, pv_lst, q1 = fn_pc(dri_lst, ri_lst, pv_lst, pi_lst, ptlc_v_lst,
                                 ptlc_s_lst, pv, divr, d0, no_surf, sign)
    while abs(q1) > 0.0000005:
        d0, q0, d1, q1 = newton_step(d0, q0, d1, q1)
        t_pv_lst, pv_lst, q1 = fn_pc(dri_lst, ri_lst, pv_lst, pi_lst,
                                     ptlc_v_lst, ptlc_s_lst, pv, divr, d0,
                                     no_surf, sign)
    return t_pv_lst, pv_lst


def fn_pc(dri_lst, ri_lst, pv_lst, pi_lst, ptlc_v_lst, ptlc_s_lst, pv, divr,
          delta, no_surf, sign):
    """
    対応するBASICコード：*PC
    """
    t_pv_lst = []
    o_pv_lst = []
    if dri_lst[0] == 0:
        t_pv_lst.append(0)
        t_pv_lst.append(pv)
    else:
        t_pv_lst.append(fn_p(divr, pv, delta))
        t_pv_lst.append(fn_p(divr, pv, delta) * divr)
    if no_surf == 4:
        o_pv_lst.extend(
            list(
                fn_ps(ptlc_s_lst[0], dri_lst[0], dri_lst[1], t_pv_lst[0],
                      pi_lst[1], ptlc_v_lst[0], pv_lst[0], pv_lst[1])))
        o_pv_lst.extend(
            list(
                fn_ps(ptlc_s_lst[1], dri_lst[2], dri_lst[3], t_pv_lst[1],
                      pi_lst[3], ptlc_v_lst[1], pv_lst[2], pv_lst[3])))
        o_q = (delta - pi_lst[1] * pv_lst[0] / t_pv_lst[0] -
               pi_lst[3] * pv_lst[3] / t_pv_lst[1] - pi_lst[2])
    elif abs(ri_lst[0]) == 1:
        o_pv_lst.append(pv_lst[0])
        o_pv_lst.extend(
            list(
                fn_ps(ptlc_s_lst[0], dri_lst[1], dri_lst[2], t_pv_lst[1],
                      pi_lst[2], ptlc_v_lst[0], pv_lst[1], pv_lst[2])))
        o_pv_lst.append(pv_lst[3])
        o_q = delta - pi_lst[2] * pv_lst[2] / t_pv_lst[1] - pi_lst[1]
    elif abs(ri_lst[1]) == 1:
        o_pv_lst.extend(
            list(
                fn_ps(ptlc_s_lst[0], dri_lst[0], dri_lst[1], t_pv_lst[0],
                      pi_lst[1], ptlc_v_lst[0], pv_lst[0], pv_lst[1])))
        o_pv_lst.append(pv_lst[2])
        o_pv_lst.append(pv_lst[3])
        o_q = delta - pi_lst[1] * pv_lst[0] / t_pv_lst[0] - pi_lst[0]
    else:
        o_pv_lst.append((t_pv_lst[0] - (sign - ri_lst[0]) * ptlc_v_lst[0]) /
                        (1 - pi_lst[1] * (sign - ri_lst[0]) * ptlc_v_lst[0]))
        o_pv_lst.append(pv_lst[1])
        o_pv_lst.append((t_pv_lst[1] - (ri_lst[1] - sign) * ptlc_v_lst[0]) /
                        (1 - pi_lst[2] * (ri_lst[1] - sign) * ptlc_v_lst[0]))
        o_pv_lst.append(pv_lst[3])
        o_q = (delta - pi_lst[1] * pv_lst[0] / t_pv_lst[0] -
               pi_lst[3] * pv_lst[3] / t_pv_lst[1] - pi_lst[2])
    return t_pv_lst, o_pv_lst, o_q


def fn_ps(ptlc_s, dri0, dri1, pv, pi, ptlc_v, o_pv0, o_pv1):
    """
    対応するBASICコード：*PS
    """
    if ptlc_s == 0:
        return dri0 * ptlc_v, (pv - o_pv0) / (1 - pi * o_pv0)
    elif ptlc_s == 1:
        return dri1 * ptlc_v, (pv - o_pv1) / (1 - pi * o_pv1)
    elif ptlc_s == 2:
        return fn_p(ptlc_v, pv, pi), fn_p(ptlc_v, pv, pi) * ptlc_v
    elif ptlc_s == 3:
        return fn_p(ptlc_v, pv, pi), fn_p(ptlc_v, pv, pi) * ptlc_v
    else:
        raise IllegalNumberError()


def fn_p(v_a, v_p, v_e):
    """
    二次方程式の解の公式の変形
    対応するBASICコード：function fnP
    """
    return 2 * v_p / (1 + v_a) / (
        1 + math.sqrt(1 - 4 * v_p * v_e * v_a / (1 + v_a)**2))


def newton_step(*pos):
    """
    ニュートン法での計算の際に用いる
    対応するBASICコード：*ND
    """
    x = (pos[0] * pos[1] - pos[2] * pos[3]) / (pos[1] - pos[3])
    return x, pos[3], pos[0], pos[3]


def rev(value):
    """
    逆数を計算する（1/0 = 0とする）
    """
    if value == 0:
        return 0
    return 1 / value


def get_stringlist(config, category, item):
    """
    configで指定されたデータファイルのitemで指定された項目が
    String型データのカンマ区切りで構成されているときに
    String型のリストとして読み込む
    """
    if config.get(category, item, fallback='') == '':
        return None
    return [
        e.strip() for e in config.get(category, item, fallback='').split(',')
    ]


def get_floatlist(config, category, item):
    """
    configで指定されたデータファイルのitemで指定された項目が
    Float型データのカンマ区切りで構成されているときに
    Float型のリストとして読み込む
    """
    if config.get(category, item, fallback='') == '':
        return None
    return [
        float(e.strip())
        for e in config.get(category, item, fallback='').split(',')
    ]


def get_intlist(config, category, item):
    """
    configで指定されたデータファイルのitemで指定された項目が
    Int型データのカンマ区切りで構成されているときに
    Int型のリストとして読み込む
    """
    if config.get(category, item, fallback='') == '':
        return None
    return [
        int(e.strip())
        for e in config.get(category, item, fallback='').split(',')
    ]


def import_glass_data(config):
    """
    configで指定されたデータファイルからガラスデータを読み込む
    対応するBASICコード：*IN4
    """
    global GLASS_DATA
    GLASS_DATA = {}
    data_index = ['dispersion', 'nd', 'nC', 'nF', 'density']
    for key, value in config['Glass Data'].items():
        prepared_value = [float(e.strip()) for e in value.split(',')]
        GLASS_DATA[key.upper()] = dict(zip(data_index, prepared_value))
