import matplotlib.pyplot as plt
import numpy as np
import csv
import math
import struct

# GPS MESSAGE [HDR|DATA|CRC]
# HDR: 3 SYNC bytes + 25 header info
# DATA: variable
# CRC: 4 bytes

##typedef struct {
##	char sync[3];
##	unsigned char headerLength;
##	unsigned short messageId;
##	char messageType;
##	unsigned char portAddress;
##	unsigned short messageLength;
##	unsigned short sequence;
##	unsigned char idleTime;
##	unsigned char timeStatus;
##	unsigned short weekNumber;
##	unsigned long gpsMsec;
##	unsigned long receiverStatus;
##	unsigned short reserved;
##	unsigned short receiverVersion;
##} GPS_Header_Typedef;

# GPS header
sync1 = 0xAA        
sync2 = 0xFF        
sync3 = 0x12        
headerLength = 28   
messageId = 241
messageType = 0         # NOT IMPLEMENTED
portAddress = 0         # NOT IMPLEMENTED
messageLength = 0       ## TO COMPLETE AT END ##
sequence = 0            # NOT IMPLEMENTED
idleTime = 0            # NOT IMPLEMENTED
timeStatus = 0          # NOT IMPLEMENTED
weekNumber = 0          # NOT IMPLEMENTED
gpsMsec = 0             # NOT IMPLEMENTED
receiverStatus = 0      # NOT IMPLEMENTED
reserved = 0            # NOT IMPLEMENTED
receiverVersion = 0     # NOT IMPLEMENTED

GPS_BESTXYZ_OFFSET_SOLSTAT = 28+0
GPS_BESTXYZ_OFFSET_POSITION = 28+8
GPS_BESTXYZ_OFFSET_POSXL = 28+8
GPS_BESTXYZ_OFFSET_POSXH = 28+12
GPS_BESTXYZ_OFFSET_POSYL = 28+16
GPS_BESTXYZ_OFFSET_POSYH = 28+20
GPS_BESTXYZ_OFFSET_POSZL = 28+24
GPS_BESTXYZ_OFFSET_POSZH = 28+28
GPS_BESTXYZ_OFFSET_POSSTDDEV = 28+32
GPS_BESTXYZ_OFFSET_PXSTD = 28+32
GPS_BESTXYZ_OFFSET_PYSTD = 28+36
GPS_BESTXYZ_OFFSET_PZSTD = 28+40
GPS_BESTXYZ_OFFSET_VELOCITY = 28+52
GPS_BESTXYZ_OFFSET_VELXL = 28+52
GPS_BESTXYZ_OFFSET_VELXH = 28+56
GPS_BESTXYZ_OFFSET_VELYL = 28+60
GPS_BESTXYZ_OFFSET_VELYH = 28+64
GPS_BESTXYZ_OFFSET_VELZL = 28+68
GPS_BESTXYZ_OFFSET_VELZH = 28+72
GPS_BESTXYZ_OFFSET_VELSTDDEV = 28+76
GPS_BESTXYZ_OFFSET_VXSTD = 28+76
GPS_BESTXYZ_OFFSET_VYSTD = 28+80
GPS_BESTXYZ_OFFSET_VZSTD = 28+84
GPS_BESTXYZ_OFFSET_NUMTRACK = 132
GPS_BESTXYZ_OFFSET_NUMUSED = 133

##typedef struct {
##	char position_solution_status[4];
##	char position_type[4];
##	double position_x;
##	double position_y;
##	double position_z;
##	float position_x_sd;
##	float position_y_sd;
##	float position_z_sd;
##	char velocity_solution_status[4];
##	char velocity_type[4];
##	double velocity_x;
##	double velocity_y;
##	double velocity_z;
##	float velocity_x_sd;
##	float velocity_y_sd;
##	float velocity_z_sd;
##      char station_id[4];
##      float v_latency;
##      float diff_age;
##      float sol_age;
##      uchar sv_tracked;
##      uchar sv_sol_used;
##      uchar gg_l1;
##      uchar soln_multi_sv;
##      char reserved;
##      char ext_sol_stat;
##      char gal_bei_sig;
##      char gps_glon_sig;
##      char[4] CRC;
##} GPS_Data_BESTXYZ;

# GPS data
position_solution_status = 0
position_type = 0
position_x = 3988310.227
position_y = 5498966.572
position_z = 900.55879
position_x_sd = 5
position_y_sd = 5
position_z_sd = 5
velocity_solution_status = 0
velocity_type = 0
velocity_x = -3290.032738
velocity_y = 2357.65282
velocity_z = 6496.623475
velocity_x_sd = 10
velocity_y_sd = 10
velocity_z_sd = 10
station_id = b'    '
v_latency = 0
diff_age = 0
sol_age = 0
sv_tracked = 10
sv_sol_used = 10
gg_l1 = 0
soln_multi_sv = 0
reserved = 0
ext_sol_stat = b' '
gal_bei_sig = b' '
gps_glon_sig = b' '
CRC = b'ABCD'

SOLSTAT = 0
NUMTRACK = 10
NUMUSED = 10


VELSDX = 10.0
VELSDY = 10.0
VELSDZ = 10.0

header = struct.pack('BBBBHBBHHBBHLLHH', sync1, sync2, sync3
            ,headerLength,messageId,messageType,portAddress
            ,messageLength,sequence,idleTime,timeStatus,weekNumber
            ,gpsMsec,receiverStatus,reserved,receiverVersion)

data = struct.pack('IIdddfffIIdddfff4sfffBBBBB1s1s1s4s'
            ,position_solution_status,position_type
            ,position_x,position_y,position_z
            ,position_x_sd,position_y_sd,position_z_sd
            ,velocity_solution_status,velocity_type
            ,velocity_x,velocity_y,velocity_z
            ,velocity_x_sd,velocity_y_sd,velocity_z_sd
            ,station_id
            ,v_latency,diff_age,sol_age
            ,sv_tracked,sv_sol_used
            ,gg_l1,soln_multi_sv
            ,reserved
            ,ext_sol_stat,gal_bei_sig,gps_glon_sig,CRC)

#struct.unpack('IIdddfffIIdddfff4sfffBBBBB1s1s1s4s',b'\x00\x00\x00\x00\x00\x00\x00\x00\x04V\x0e\x1d\xabmNA\xe3\xa5\x9b\xa4\x15\xfaTA\xa8:\xe4fx$\x8c@\x00\x00\xa0@\x00\x00\xa0@\x00\x00\xa0@\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\xac\xfe\x08\xc3\x10\xb4\xa9\xc0YLl>Nk\xa2@\xe0\xbe\x0e\x9c\x9f`\xb9@\x00\x00 A\x00\x00 A\x00\x00 A    \x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\n\n\x00\x00\x00   ABCD')
