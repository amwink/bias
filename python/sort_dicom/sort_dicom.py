# to run from the python prompt:
# >>> exec(open('dicoms_prot.py').read());

import pydicom;
import json;

protfile="dicoms_prot.json";
dcmslist="dicoms_list.txt";

debug_level = 0;

# load the protocol definitions
with open ( protfile, "r" ) as inputfile:
    protocol_data = json.load ( inputfile );

# make access to the sets easier:
protocol_sets = { s['Name']: s['Values'] for s in protocol_data [ 'sets' ] };
    
# load the protocol definitions
counter = 0;
with open ( dcmslist, "r" ) as inputfile:

    dcmnames = inputfile.readlines();

    for filename in dcmnames:

        if debug_level > 0:
            print ( filename )
            
        protfound = "unknown";
        
        # read dicom record (header of 1st slice in series)
        dcmrecord = pydicom.dcmread ( filename.strip() )

        # iterate protocol list until the correct protocol
        # (one whose parameter values all match) is found
        for p in protocol_data["protocols"]:

            tested_keys  = 0;
            correct_keys = 0;
                    
            for k in p.keys():
                
                if not ( k == "Name" ):
                
                    if debug_level > 0:
                        print ( 'protocol {}: key {}'.format( p["Name"]["free"], k ) );                
                    
                    if "set" in p[k] and (correct_keys == tested_keys):     # if this parameter is of the "set" type
                        tested_keys += 1;                                   #
                        if k in dcmrecord:                                  # if the dicoom file has this parameter and
                            if debug_level > 1:                             #
                                print ( 'value: {}'.format ( dcmrecord[k].value ) );
                            if ( p[k]["set"] in protocol_sets [k] ):        # if the value is part of the valid set of values
                                if ( p[k]["set"] in dcmrecord[k].value ):   # if the value is the one for this protocol
                                    correct_keys += 1;                      # increase number of fits

                        if debug_level > 1:
                            if tested_keys == correct_keys:
                                print ( 'set: value {} matches for {}'.format ( dcmrecord[k].value, p["Name"]["free"] ) );
                            else:
                                print ( 'set: value for {} not matched ({}) in set {}'.format ( k, p[k]["set"], protocol_sets[k] ) );
                                break;
                        
                    if "min" in p[k] and (correct_keys == tested_keys):                   # if this parameter is a lower bound
                        tested_keys += 1;                                                 #
                        if k in dcmrecord:                                                # see if it is in the dicom file
                            if debug_level > 1:
                                print ( 'value: {}'.format ( dcmrecord[k].value ) );
                            if ( float ( dcmrecord[k].value ) >= float ( p[k]["min"] ) ): # check if the value in the dicom is higher
                                correct_keys += 1;                                        # increase number of fits

                        if debug_level > 1:
                            if tested_keys == correct_keys:
                                print ( 'min: value {} matches for {}'.format ( dcmrecord[k].value, p["Name"]["free"] ) );
                            else:
                                print ( 'min: value for {} not found or too low (min {})'.format ( k , p[k]["min"] ) );
                                break;
                            
                    if "max"  in p[k] and (correct_keys == tested_keys):                  # if this parameter is a lower bound
                        tested_keys += 1;                                                 #
                        if k in dcmrecord:                                                # see if it is in the dicom file
                            if debug_level > 1:
                                print ( 'value: {}'.format ( dcmrecord[k].value ) );
                            if ( float ( dcmrecord[k].value ) <= float ( p[k]["max"] ) ): # check if the value in the dicom is lower
                                correct_keys += 1;                                        # increase number of fits

                        if debug_level > 1:
                            if correct_keys == tested_keys:
                                print ( 'max: value {} matches for {}'.format ( dcmrecord[k].value, p["Name"]["free"] ) );
                            else:
                                print ( 'max: value for {} not found or too high (max {})'.format ( k, p[k]["max"] ) );
                                
                     # end if k not "Name"
                   
                if debug_level >1 and correct_keys < tested_keys:
                    print ( "tested: {}, correct: {}, continuing ...".format ( tested_keys, correct_keys ) );
                    break;
                     
                # end loop over k

            if debug_level > 0:
                print ( "protocol {}, fitting keys: {}".format ( p, correct_keys ) );
            if ( correct_keys >= len ( p.keys() ) -2 ):
                protfound = p["Name"]["free"];
                break;
                        
            # end loop over p

        # prepare for next file
        if not ( protfound == "unknown" ):
            print ( '{}: {}'.format ( protfound, filename.strip() ) );
        counter += 1;

        if debug_level > 0:
            try:
                input("Press enter to continue")
            except SyntaxError:
                pass

        # end for filenames

    # end with inputfile

    
