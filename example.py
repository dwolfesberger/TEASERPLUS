from teaser.project import Project
import teaser.data.input.citygml_input as cg
import teaser.data.output.citygml_output as cg_out
import teaser.logic.utilities as utilities
import os



def example_citygml():

    prj1 = Project(load_data=True)
    prj1.name = "Teaserworkshop_1"
con
    """TEASER - Archetypes"""
    # prj1.add_non_residential(method='bmvbs', name='Office', year_of_construction=1980, net_leased_area=240,
    #                      usage='office', number_of_floors=2, height_of_floors=3.2588)
    #
    # prj1.add_residential(method='iwu', name='SFD', year_of_construction=1980, net_leased_area=240,
    #                      usage='single_family_dwelling', number_of_floors=2, height_of_floors=3.2588)
    #
    # prj1.add_residential(method='tabula_de', name='SFH', year_of_construction=1980, net_leased_area=240,
    #                      usage='single_family_house', number_of_floors=2, height_of_floors=3.2588)

    """CityGML LoD3 Model """
    cg.load_gml('../../gmlfiles/FZK-House-LoD0-LoD4/FZK-Haus-LoD3-KIT-IAI-KHH-B36-V1.gml', prj1, method='iwu')

    """Multiple CityGML Residential Buildings / BuildingParts"""
    # gml_copy_list_LoD2 = cg.choose_gml('../teaserplus/gmlfiles/LoD2_Aachen_UncertaintyProject/Aachen_LoD2_31001_1010.gml', prj1,
    #                                    bldg_names=['DENW39AL1000MjqW',
    #                                                'DENW39AL1000MjeA',
    #                                                'DENW39AL1000Mje8',
    #                                                'DENW39AL10005X5J',
    #                                                'DENW39AL10005VOB'])
    # cg.load_gml('../teaserplus/gmlfiles/LoD2_Aachen_UncertaintyProject/Aachen_LoD2_31001_1010.gml', prj1, method='iwu',
    #             chosen_gmls=gml_copy_list_LoD2)

    """EnergyADE-files"""
    # cg.load_gmlade('../teaserplus/gmlfiles/FZKHouseLoD2-ADE.gml', prj1)

    """Additional Parameters"""
    path = utilities.get_default_path()
    prj1.number_of_elements_calc = 4
    prj1.calc_all_buildings(raise_errors=False)
    prj1.export_parameters_txt()
    prj1.weather_file_path = utilities.get_full_path(
        os.path.join(
            "data",
            "input",
            "inputdata",
            "weatherdata",
            "DRYCOLDTMY.mos"))
    prj1.export_aixlib()

    """GML out"""
    # prj1.save_citygml(path="C:/Users/You/teaserplus/gmlfiles/ADE out/TEASER_workshop")

    # prj1.save_citygml(path="C:/Users/You/teaserplus/gmlfiles/ADE out/Aachen_LoD2_31001_1010_ADE"
    #                   , gml_copy=gml_copy_list_LoD2)


    '''Start the simulation and save the results'''
    # s.simulate(path=path, prj=prj1, loading_time=None,
    #            result_path=f'C:/Users/MaxPaine33/PycharmProjects/teaserplus/gmlfiles/'
    #                        f'ADE out/TEASER_workshop/{prj1.name}_results.csv')


if __name__ == '__main__':
    example_citygml()
    print("ShouldWork!")