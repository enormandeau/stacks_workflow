#rename 's/fasta/fa/' *.fasta
for i in `ls -1 *.unwrapped`; do perl -ne 'if (m/>/) {print} else {print(substr $_, 12)}' $i > $i".trimmed"; done &
#rename 's/\.fa/.fa.bkp/' *.fa
#rename 's/\.fa\.trimmed/.fa/' *.trimmed

