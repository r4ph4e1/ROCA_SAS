#!/usr/bin/bash

PUBKEY="key.pub"
KEYSIZE=512
FILE="secret.txt"
OUT="secret.enc"
PRIVKEY="key.priv"

echo "Dieser Text ist streng geheim" > $FILE

echo "Generate Key"
python2 keygen.py $PUBKEY $KEYSIZE

echo "Generating Key-File"
dd if=/dev/urandom of=key.bin bs=35 count=1

# Verschlüsseln eines Files mit AES:
echo "Encrypting the file $FILE with AES-256-CBC and the Key-File as Password"
openssl enc -aes-256-cbc -salt -in $FILE -out  $OUT -pass file:./key.bin

# Verschlüsseln des Key-Files "key.bin" mit unserem RSA Key
echo "Encrypting the Key-File with the RSA-Public-Key: $PUBKEY"
openssl rsautl -encrypt -inkey $PUBKEY -pubin -in key.bin -out key.bin.enc

# Löschen des unverschlüsselten Dokuments
echo "Deleting unencrypted Files: $FILE and key.bin"
rm $FILE key.bin

sage roca.sage $PUBKEY $PRIVKEY

#Entlüsseln des Key-Files "key.bin" mit dem Faktorisieren Private Key
echo "Decrypting the Key-File with the factorized RSA-Private-Key"
openssl rsautl -decrypt -inkey $PRIVKEY -in key.bin.enc -out key.bin

# Entschlüsseln der AES verschlüsselten Datei mithilfe des entschlüsselten Key-Files:
echo "Decrypting the File $FILE with the decrypted Key-File"
openssl enc -d -aes-256-cbc -in $FILE -out $OUT -pass file:./key.bin

# Löschen der verschlüsselten Dokumente
echo "Deleting Encrypted Files: $FILE and key.bin.enc"
rm $FILE "key.bin.enc"


