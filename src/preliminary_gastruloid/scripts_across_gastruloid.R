assign_type_to_number = function(type){

  if (type == "Prog"){
    return(0)
  }else if (type == " SMD"){
    return(1)
  }else if (type == " NMP"){
    return(2)
  }else if (type == " PhMD"){
    return(3)
  }else if(type == " PxMD"){
    return(4)
  }else if(type == " SC"){
    return(5)
  }else{
    stop("Type not found!")
  }
}
