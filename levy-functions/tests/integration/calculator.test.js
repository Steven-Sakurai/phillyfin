const handler=require('../../lambda/calculatorLambda')
const createEvent=(data, parameters)=>({
    body:JSON.stringify(data),
    pathParameters:parameters
})
it('correctly returns heston price', (done)=>{
    //http://ta.twi.tudelft.nl/mf/users/oosterle/oosterlee/COS.pdf pg 15
    const r=0.0
    const T=1.0
    const S0=100
    const b=.0398
    const a=1.5768
    const c=.5751
    const rho=-.5711
    const v0=.0175

    //convert parameters
    const sigma=Math.sqrt(b)
    const speed=a
    const v0Hat=v0/b
    const adaV=c/Math.sqrt(b)

    const parameters={
        numU:8,
        r,
        T,
        S0,
        sigma, 
        lambda:0,
        muJ:0,
        sigJ:0,
        speed,
        v0:v0Hat,
        adaV,
        rho,
        k:[100]
    }
    const event=createEvent(parameters, {
        optionType:'call',
        sensitivity:'price',
        algorithm:'fangoost'
    })
    return handler.calculator(event, {}, (err, val)=>{
        console.log(val.body)
        const parsedVal=JSON.parse(val.body)
        expect(parsedVal[1].value).toBeCloseTo(5.78515545, 3)
        done()
    })
})
it('correctly returns merton price', (done)=>{
    //https://www.upo.es/personal/jfernav/papers/Jumps_JOD_.pdf pg 8
    const r=.1
    const T=.5
    const S0=38
    const sigJ=Math.sqrt(.05)
    const sigma=sigJ
    const muJ=-sigJ*sigJ*.5
    const k=35
    const lambda=1
    
    const parameters={
        numU:8,
        r,
        T,
        S0,
        sigma, 
        lambda,
        muJ,
        sigJ,
        speed:0,
        v0:1,
        adaV:0,
        rho:0,
        k:[k]
    }
    const event=createEvent(parameters, {
        optionType:'call',
        sensitivity:'price',
        algorithm:'fangoost'
    })
    return handler.calculator(event, {}, (err, val)=>{
        console.log(val.body)
        const parsedVal=JSON.parse(val.body)
        expect(parsedVal[1].value).toBeCloseTo(5.9713, 3)
        done()
    })
})
it('correctly returns VaR', (done)=>{
    //https://github.com/phillyfan1138/levy-functions/issues/27
    const r=.004
    const T=.187689
    const S0=191.96
    const sigJ=.220094
    const sigma=.3183
    const muJ=-.302967
    const lambda=.204516
    const speed=2.6726
    const v0=.237187
    const rho=-.182754
    const adaV=0
    
    const parameters={
        numU:8,
        r,
        T,
        S0,
        sigma, 
        lambda,
        muJ,
        sigJ,
        speed,
        v0,
        adaV,
        rho,
        quantile:.01
    }
    const event=createEvent(parameters, {
        densityType:'var'
    })
    return handler.density(event, {}, (err, val)=>{
        console.log(val.body)
        const parsedVal=JSON.parse(val.body)
        expect(parsedVal.VaR).toBeCloseTo(.261503, 3)
        done()
    })
})
