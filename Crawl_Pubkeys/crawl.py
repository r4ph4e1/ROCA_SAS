

import subprocess
import netaddr
import ipaddress
import datetime

'''
TODO
- Ip-Ranges Parsen
- Ranges nach offenen Ports absuchen, ergebnisse je range in eigenes File schreiben
- IP's in file speichern
- entwedet nmap keyscan oder ssh-keyscan durchfuehren
- Informationen speichern (IP, Pubkey auth verfuegbar Ja/Nein)
- keys als pem speichern
- mit Fingerprinting Tool ueberpruefen
- Key Kategorisieren nach keylaengen etc.
- Gewonnene Informationen mergen (Offene IPs, anfaellig, nicht anfaellig, keylaengen)
'''

def nmap_scan(ip_range):
    print('NMAP: Scanning Range: ' + ip_range)
    #nmap -n -sn 192.0.2.0/24 -oG - | awk '/Up$/{print $2}'
    args = "nmap -sV " + ip_range + " -p 22 -oG -"
    args = args.split(" ")
    p = subprocess.Popen(args, stdout=subprocess.PIPE)
    output = subprocess.check_output(('awk', '/Up$/{print $2}'), stdin=p.stdout)
    print(p.stdout)
    p.wait()

    output = str(output, 'utf-8').rstrip()
    ips = output.split('\n')

    f = open('output/' + 'IPs_' + ip_range.split('/', 1)[0] + '.txt', 'w+')
    f.write(output)
    f.close()

    get_pubkeys('output/' + 'IPs_' + ip_range.split('/', 1)[0] + '.txt')
    return ips

def get_pubkeys(ip_file):

    print('Searching PubKeys for: ' + ip_file)

    ips = []

    with open(ip_file, 'r') as fp:
        for line in fp:
            line = line.split('\n', 1)[0]
            ips.append(line)

    print(ips)

    for ip in ips:
        args = "ssh-keyscan -t rsa " + ip
        args = args.split(" ")

        p = subprocess.Popen(args, stdout=subprocess.PIPE)
        output = p.communicate()[0]
        output = str(output, 'utf-8').rstrip()
        output = output.split(' ', 1)[-1]

        if output:
            f = open('output/' + ip + '.pub', 'w+')
            f.write(output)
            f.close()


def parse_ranges(range_file):
    ranges = []
    ctr = 0
    with open(range_file, 'r') as fp:
        for line in fp:
            line = line.split('\n', 1)[0]
            start, end = line.split(',')
            #print(start, end)
            cidr = netaddr.iprange_to_cidrs(start, end)

            if len(cidr) > 1:
                for i in range(0, len(cidr)):
                    ranges.append(cidr[i])
                    range_size = netaddr.IPNetwork(cidr[i])
            else:
                ip_range = ''.join(str(c) for c in cidr)
                ranges.append(ip_range)

                range_size = netaddr.IPNetwork(ip_range)

            ctr += range_size.size # Anzahl an IP Adressen aller Ranges

    f = open('ranges.txt', 'w+')
    for ip_range in ranges:
        f.write(str(ip_range) + '\n')
    f.close()


#get_pubkeys('IPs_217.199.80.0.txt')
#parse_ranges('start_end_ips.txt')
num_ranges = sum(1 for line in open('ranges.txt'))
count = 1
timestamp = datetime.datetime.now()

f = open('ranges.txt', 'r')

for line in f:
    print('Start Time: %s' % (str(timestamp)))
    line = line.split('\n', 1)[0]
    print("Scanning Range %d of %d" % (count, num_ranges))
    nmap_scan(line)
    count += 1
    used_time = datetime.datetime.now() - timestamp
    print('Used Time: %s' % (str(used_time)))
    timestamp = datetime.datetime.now()
