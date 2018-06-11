#!/usr/bin/bash

while getopts p:f:o: option
do
	case "${option}" in
		f) FILE=${OPTARG};;
		p) PUBKEY=${OPTARG};;
		o) OUT=${OPTARG};;
	esac
done

# Generieren eines passwort-files der Größe 35 (geht nicht viel größer da der 512-bit Key zu klein ist und er einen größeren nicht verschlüsseln kann)
echo "Generating Key-File"
openssl rand -base64 35 -out key.bin


# Verschlüsseln eines Files mit AES:
echo "Encrypting the file $FILE with AES-256-CBC and the Key-File as Password"
openssl enc -aes-256-cbc -salt -in $FILE -out  $OUT -pass file:./key.bin

# Verschlüsseln des Key-Files "key.bin" mit unserem RSA Key
echo "Encrypting the Key-File with the RSA-Public-Key: $PUBKEY"
openssl rsautl -encrypt -inkey $PUBKEY -pubin -in key.bin -out key.bin.enc

# Löschen des unverschlüsselten Dokuments
echo "Deleting unencrypted Files: $FILE and key.bin"
rm $FILE key.bin
