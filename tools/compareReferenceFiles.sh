if [ -f u_ref.dat ]
then
  echo "number of lines of u.dat $(cat u.dat | wc -l) | $(cat u_ref.dat | wc -l)"
  #diff u.dat u_ref.dat
  counter=$(diff -y --suppress-common-lines u.dat u_ref.dat | wc -l)
  if [[ $counter > 0 ]]
  then
    echo "u.dat has $counter differences"
  fi
fi
if [ -f v_ref.dat ]
then
  echo "number of lines of v.dat $(cat v.dat | wc -l) | $(cat v_ref.dat | wc -l)"
  #diff v.dat v_ref.dat
  counter=$(diff -y --suppress-common-lines v.dat v_ref.dat | wc -l)
  if [[ $counter > 0 ]]
  then
    echo "v.dat has $counter differences"
  fi
fi
if [ -f w_ref.dat ]
then
  echo "number of lines of w.dat $(cat w.dat | wc -l) | $(cat w_ref.dat | wc -l)"
  #diff w.dat w_ref.dat
  counter=$(diff -y --suppress-common-lines w.dat w_ref.dat | wc -l)
  if [[ $counter > 0 ]]
  then
    echo "w.dat has $counter differences"
  fi
fi
if [ -f p_ref.dat ]
then
  echo "number of lines of p.dat $(cat p.dat | wc -l) | $(cat p_ref.dat | wc -l )"
  #diff p.dat p_ref.dat
  counter=$(diff -y --suppress-common-lines p.dat p_ref.dat | wc -l)
  if [[ $counter > 0 ]]
  then
    echo "p.dat has $counter differences"
  fi
fi
if [ -f T_ref.dat ]
then
  echo "number of lines of T.dat $(cat T.dat | wc -l) | $(cat T_ref.dat | wc -l)"
  #diff T.dat T_ref.dat
  counter=$(diff -y --suppress-common-lines T.dat T_ref.dat | wc -l)
  if [[ $counter > 0 ]]
  then
    echo "T.dat has $counter differences"
  fi
fi

