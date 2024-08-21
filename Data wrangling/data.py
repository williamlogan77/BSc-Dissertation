from . import extras, x_time_unpacker, mass_unpacker, flux_unpacker


def get_prepared_data(x_time_path='/home/william/Desktop/Uni/ppn/nuppn/frames/ppn/i_process/x-time.dat', ppn_path='../Data_files/iniab1.4E-02As09.ppn', flux_path='../Data_files/flux_00010.DAT'):
    with open(x_time_path, 'r') as data_file:
        lines = data_file.readlines()

    with open(ppn_path, 'r') as data_file:
        mass_lines = data_file.readlines()

    with open(flux_path, 'r') as data_file:
        flux_lines = data_file.readlines()

    data_dict = x_time_unpacker.get_data_from_lines(lines)
    data_dict = mass_unpacker.unpack_mass(mass_lines, data_dict)
    # data_dict = flux_unpacker.unpack_flux(flux_lines, data_dict)
    return extras.fill_blanks(data_dict)
